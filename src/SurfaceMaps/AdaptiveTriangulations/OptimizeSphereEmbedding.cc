/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Dörte Pieper, Patrick Schmidt
 */

#include "OptimizeSphereEmbedding.hh"

#include <TinyAD/Support/OpenMesh.hh>
#include <TinyAD/ScalarFunction.hh>
#include <TinyAD/Utils/NewtonDirection.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>
#include <TinyAD/Utils/LineSearch.hh>
#include <TinyAD/Utils/Timer.hh>
#include <TinyAD/Utils/Out.hh>
#include <SurfaceMaps/AdaptiveTriangulations/Helpers.hh>
#include <SurfaceMaps/Misc/DistortionEnergy.hh>
#include <SurfaceMaps/Misc/ConstantCurvatureEmbedding.hh>

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
        SVH v,
        const ExternalProperty<VH, Vec3d>& _embedding,
        const ExternalProperty<VH, Vec3d>& _B1,
        const ExternalProperty<VH, Vec3d>& _B2)
{
    // Evaluate target point in 3D ambient space and project to sphere via normalization.
    return (_embedding[v]  + v_tang[0] * _B1[v] + v_tang[1] * _B2[v]).normalized();
}

/**
 * Signed area of spherical triangle
 */
double signed_spherical_area(
        const Vec3d& _a,
        const Vec3d& _b,
        const Vec3d& _c)
{
    // Plane normals
    Vec3d n_ab = _a.cross(_b).normalized();
    Vec3d n_bc = _b.cross(_c).normalized();
    Vec3d n_ca = _c.cross(_a).normalized();

    // Inner angles
    double alpha = M_PI - acos(n_ab.dot(n_bc));
    double beta = M_PI - acos(n_bc.dot(n_ca));
    double gamma = M_PI - acos(n_ca.dot(n_ab));

    // Check orientation
    double sign = flipped_or_degenerate(_a, _b, _c, Spherical) ? -1.0 : 1.0;

    // Signed spherical area
    return sign * (alpha + beta + gamma - M_PI);
}

double map_degree(
        const TriMesh& _mesh,
        const ExternalProperty<VH, Vec3d>& _embedding)
{
    double total_area = 0.0;
    for (auto f : _mesh.faces())
    {
        SVH va, vb, vc;
        handles(_mesh, f, va, vb, vc);
        Eigen::Vector3d a = _embedding[va];
        Eigen::Vector3d b = _embedding[vb];
        Eigen::Vector3d c = _embedding[vc];
        total_area += signed_spherical_area(a, b, c);
    }

    return total_area / 4.0 / M_PI;
}

double map_degree(
        const TriMesh& _mesh,
        const VecXd& _x,
        const ExternalProperty<VH, Vec3d>& _embedding,
        const ExternalProperty<VH, Vec3d>& _B1,
        const ExternalProperty<VH, Vec3d>& _B2)
{
    ExternalProperty<VH, Vec3d> tmp_embedding(_mesh);
    for (auto v : _mesh.vertices())
        tmp_embedding[v] = retract(Vec2d(_x.segment(2 * v.idx(), 2)), v, _embedding, _B1, _B2);

    return map_degree(_mesh, tmp_embedding);
}

}

bool sphere_embedding_bijective(
        const TriMesh& _mesh,
        const ExternalProperty<VH, Vec3d>& _embedding)
{
    if (count_flipped(_mesh, _embedding, Spherical) != 0)
        return false;

    if (round(map_degree(_mesh, _embedding)) != 1)
        return false;

    return true;
}

void optimize_sphere_embedding(
        const TriMesh& _mesh,
        ExternalProperty<VH, Vec3d>& _embedding,
        const int _n_iters,
        const OptimizeSphereEmbeddingSettings& _settings,
        std::function<void()> _callback)
{
    const double surface_area = total_area(_mesh);

    // Compute an orthonormal tangent-space basis at each vertex.
    ExternalProperty<VH, Vec3d> B1(_mesh);
    ExternalProperty<VH, Vec3d> B2(_mesh);
    compute_local_bases(_embedding, B1, B2, _mesh);

    // Set up function with a 2D tangent vector per vertex as variables.
    auto func = TinyAD::scalar_function<2>(_mesh.vertices());

    // Add triangle elements, each connecting three vertices.
    func.add_elements<3>(_mesh.faces(), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
    {
        // Evaluate element using either double or TinyAD::Double
        using T = TINYAD_SCALAR_TYPE(element);

        // Get triangle and vertex handles
        OpenMesh::SmartFaceHandle f = element.handle;
        auto va = f.halfedge().to();
        auto vb = f.halfedge().next().to();
        auto vc = f.halfedge().from();

        // Get vertex positions, each in their own 2D tangent space.
        // These are (0, 0) when computing derivatives and != (0, 0) in line search.
        Eigen::Vector2<T> a_tang = element.variables(va);
        Eigen::Vector2<T> b_tang = element.variables(vb);
        Eigen::Vector2<T> c_tang = element.variables(vc);

        // Apply retraction operator: Translate tangent vectors into 3D points on the sphere.
        Eigen::Vector3<T> a_mani = retract(a_tang, va, _embedding, B1, B2);
        Eigen::Vector3<T> b_mani = retract(b_tang, vb, _embedding, B1, B2);
        Eigen::Vector3<T> c_mani = retract(c_tang, vc, _embedding, B1, B2);

        return sphere_distortion_energy(
                    f, a_mani, b_mani, c_mani, _mesh, surface_area,
                    _settings.w_barrier, _settings.w_angle, _settings.w_area, _settings.max_edge_length_degrees);
    });

    // Variable vector: 2D tangent vector per vertex.
    // Initially, each vertex sits at the tangent-space origin.
    Eigen::VectorXd x = Eigen::VectorXd::Zero(2 * _mesh.n_vertices());

    // Optimize via Projected-Newton
    double f;
    Eigen::VectorXd g;
    Eigen::SparseMatrix<double> H_proj;
    Eigen::VectorXd d;
    TinyAD::LinearSolver solver;
    int iter = 0;
    for (; iter < _n_iters; ++iter)
    {
        // callback
        _callback();

        // Compute derivatives
        func.eval_with_hessian_proj(x, f, g, H_proj, _settings.hessian_proj_eps);
//        TINYAD_DEBUG_OUT("Energy in iteration " << iter << ": " << f);

        // Compute update direction
        d = TinyAD::newton_direction(g, H_proj, solver, _settings.w_identity);

        // Converged?
//        TINYAD_DEBUG_VAR(TinyAD::newton_decrement(d, g));
        if (TinyAD::newton_decrement(d, g) < _settings.convergence_eps)
        {
            ISM_HIGHLIGHT("Sphere optimization converged in iteration " << iter);
            break;
        }

        // Line search
        x = TinyAD::line_search(x, d, f, g, [&] (const VecXd& _x_ls)
        {
            if (round(map_degree(_mesh, _x_ls, _embedding, B1, B2)) != 1)
                return INF_DOUBLE;

            return func.eval(_x_ls);
        });

        // Re-center local bases
        for (auto v_idx : _mesh.vertices())
            _embedding[v_idx] = retract(Vec2d(x.segment(2 * v_idx.idx(), 2)), v_idx, _embedding, B1, B2);
        compute_local_bases(_embedding, B1, B2, _mesh);
        x = Eigen::VectorXd::Zero(2 * _mesh.n_vertices());
    }

    _callback();

//    TINYAD_DEBUG_OUT("Final energy: " << func.eval(x));
}

}
