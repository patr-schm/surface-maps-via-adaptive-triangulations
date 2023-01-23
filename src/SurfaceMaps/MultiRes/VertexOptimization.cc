/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Joe Jakobi, Patrick Schmidt
 */

#include "VertexOptimization.hh"

#include <SurfaceMaps/Types.hh>
#include <SurfaceMaps/Utils/Helpers.hh>
#include <SurfaceMaps/AdaptiveTriangulations/Helpers.hh>
#include <SurfaceMaps/Misc/ConstantCurvatureGeometry.hh>
#include <SurfaceMaps/Misc/ConstantCurvatureEmbedding.hh>
#include <SurfaceMaps/Misc/DistortionEnergy.hh>
#include <TinyAD/Utils/Timer.hh>
#include <TinyAD/Scalar.hh>

namespace SurfaceMaps
{

namespace
{

/**
 * Retraction operator: map from a local tangent space to the sphere.
 */
template <typename T>
Vec3<T> retract(
        const Vec2<T>& v_tang,
        const Vec3d& _p,
        const Vec3d& _b1,
        const Vec3d& _b2)
{
    // Evaluate target point in 3D ambient space and project to sphere via normalization.
    return (_p  + v_tang[0] * _b1 + v_tang[1] * _b2).normalized();
}

template <typename T>
T triangle_energy(
        const Vec2<T>& _x,
        const FH _fh,
        const VH _vh_a, // Variable
        const VH _vh_b, // Constant
        const VH _vh_c, // Constant
        const Vec3d& _b1,
        const Vec3d& _b2,
        const MultiResSphereEmbedding& _data,
        const MultiResSphereEmbeddingSettings& _settings)
{
    const Vec3d p_a_orig = _data.mesh.property(_data.ph_embedding, _vh_a);

    const Vec3<T> p_a = retract(_x, p_a_orig, _b1, _b2);
    const Vec3<double> p_b = _data.mesh.property(_data.ph_embedding, _vh_b);
    const Vec3<double> p_c = _data.mesh.property(_data.ph_embedding, _vh_c);

    return sphere_distortion_energy<T>(
                _fh, p_a, p_b, p_c,
                _data.mesh, _data.mesh_area,
                _settings.refinement_opt_settings.w_barrier,
                _settings.refinement_opt_settings.w_angle,
                _settings.refinement_opt_settings.w_area,
                _settings.refinement_opt_settings.max_edge_length_degrees);
}

double eval_one_ring(
        const Vec2d& _x,
        const VH _vh,
        const Vec3d& _b1,
        const Vec3d& _b2,
        const MultiResSphereEmbedding& _data,
        const MultiResSphereEmbeddingSettings& _settings,
        Vec2d* _g = nullptr,
        Mat2d* _H = nullptr)
{
    ISM_ASSERT_FINITE_MAT(_x);
    ISM_ASSERT_EQ(_g == nullptr, _H == nullptr);

    using ADouble = TinyAD::Double<2>;

    double f = 0.0;
    const Vec2<ADouble> x_a = ADouble::make_active(_x);
    for (auto heh : _data.mesh.voh_range(_vh))
    {
        if (heh.is_boundary())
            continue;

        VH vh_a = heh.from();
        VH vh_b = heh.to();
        VH vh_c = heh.next().to();

        if (_g)
        {
            ADouble f_a = triangle_energy<ADouble>(x_a, heh.face(), vh_a, vh_b, vh_c, _b1, _b2, _data, _settings);

            f += f_a.val;
            *_g += f_a.grad;
            *_H += f_a.Hess;
        }
        else
        {
            f += triangle_energy<double>(_x, heh.face(), vh_a, vh_b, vh_c, _b1, _b2, _data, _settings);
        }
    }

    return f;
}

double max_step(
        const TriMesh& _mesh,
        const VH _vh,
        const OpenMesh::VPropHandleT<Vec3d>& _ph_param,
        const Vec2d& _d)
{
    const Vec3d p = _mesh.property(_ph_param, _vh);
    double max_tangent_length = 0.0;
    for (auto v_neigh : _mesh.vv_range(_vh))
    {
        const Vec3d p_neigh = _mesh.property(_ph_param, v_neigh);
        const double dot = p.dot(p_neigh);

        if (dot <= 0.0)
            continue;

        const double tangent_length_to_neigh = tan(acos(dot));
        max_tangent_length = std::max(max_tangent_length, tangent_length_to_neigh);
    }

    if (max_tangent_length <= 0.0)
        max_tangent_length = tan(75.0 / 180.0 * M_PI);

    return max_tangent_length / _d.norm();
}

namespace
{

bool armijo_condition(
        const double _f_curr,
        const double _f_new,
        const double _s,
        const VecXd& _d,
        const VecXd& _g,
        const double _armijo_const)
{
    return _f_new <= _f_curr + _armijo_const * _s * _d.dot(_g);
}

// TODO: Use TinyAD::line_search
bool line_search(
        const VecXd& _x0,
        const VecXd& _g,
        const VecXd& _d,
        const double _s_max,
        const double _f,
        VecXd& _x_new,
        double& _f_new,
        std::function<double(
            const VecXd& _x)> _eval_f,
        const bool _silent)
{
    constexpr int line_search_max_iters = 64;
    constexpr double line_search_shrink = 0.8;
    constexpr double armijo_const = 1e-4; // The larger the tighter (i.e. smaller steps)

    // Also try a step size of 1.0 (if valid)
    const bool try_one = _s_max > 1.0;

    TinyAD::Timer timer("Line search", _silent);

    ISM_ASSERT_EQ(_x0.size(), _g.size());
    if (_s_max <= 0.0)
        ISM_ERROR_throw("Max step size not positive.");

    double s = _s_max;
    for (int i = 0; i < line_search_max_iters; ++i)
    {
        _x_new = _x0 + s * _d;
        _f_new = _eval_f(_x_new);
        ISM_ASSERT_NOT_NAN(_f_new);
        if (armijo_condition(_f, _f_new, s, _d, _g, armijo_const))
            return true;

        if (try_one && s > 1.0 && s * line_search_shrink < 1.0)
            s = 1.0;
        else
            s *= line_search_shrink;
    }

    if (!_silent)
        ISM_WARNING("Line search failed.");

    return false;
}

}

bool perform_step(
        const MultiResSphereEmbedding& _data,
        const VH _vh,
        const Vec3d& _b1,
        const Vec3d& _b2,
        const MultiResSphereEmbeddingSettings& _settings,
        const Vec2d& _x,
        const Vec2d& _g,
        const Vec2d& _d,
        const double& _f,
        VecXd& _x_new,
        double& _f_new)
{
    const double s_max = max_step(_data.mesh, _vh, _data.ph_embedding, _d);

    const bool success = line_search(_x, _g, _d, s_max, _f, _x_new, _f_new, [&] (const Vec2d& _x)
    {
        return eval_one_ring(_x, _vh, _b1, _b2, _data, _settings);
    }, true);
    return success;
}

}

void optimize_vertex(
        MultiResSphereEmbedding& _data,
        const VH _vh,
        const uint _max_iters,
        const MultiResSphereEmbeddingSettings& _settings)
{
    const HEH heh = _data.mesh.halfedge_handle(_vh);
    ISM_ASSERT(_data.mesh.from_vertex_handle(heh) == _vh);

    for (uint i = 0; i < _max_iters; ++i)
    {
        // Compute tangent space basis
        const Vec3d p_orig = _data.mesh.property(_data.ph_embedding, _vh);
        const Vec3d b1 = any_orthogonal(p_orig);
        const Vec3d b2 = p_orig.cross(b1);

        Vec2d x = Vec2d(0.0, 0.0); // (0, 0) in tangent-space
        Vec2d g = Vec2d::Zero(); // Init important!
        Mat2d H = Mat2d::Zero(); // Init important!

        double f = eval_one_ring(x, _vh, b1, b2, _data, _settings, &g, &H);

        ISM_ASSERT_FINITE(f);

        // Converged?
        if (g.norm() < _settings.optimization_gradient_thresh)
            break;

        // Compute update direction
        H += _settings.w_identity * Mat2d::Identity();

        const Vec2d d = -H.inverse() * g;

        // Perform Newton step or fall back to gradient descent
        VecXd x_new;
        double f_new;

        bool success = perform_step(_data, _vh, b1, b2, _settings, x, g, d, f, x_new, f_new);
        if (!success)
        {
//            ISM_WARNING("Falling back to gradient descent because line search failed.");
            success = perform_step(_data, _vh, b1, b2, _settings, x, g, -g, f, x_new, f_new);
        }

        if (!success)
        {
            ISM_WARNING("Optimization of vertex " << _vh << " in iteration " << i << " failed.");
            break;
        }

        ISM_EXPECT_LEQ(f_new, f);
        x = x_new;

        _data.mesh.property(_data.ph_embedding, _vh) = retract<double>(x, p_orig, b1, b2);

        // Converged
        if (f - f_new < _settings.optimization_improvment_thresh)
            break;
    }
}

}
