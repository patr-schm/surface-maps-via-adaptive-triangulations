/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Dörte Pieper, Patrick Schmidt
 */

#include "OptimizeMapEnergies.hh"

#include <TinyAD/Scalar.hh>
#include <TinyAD/Utils/Helpers.hh>

#include <SurfaceMaps/AdaptiveTriangulations/AdaptiveTargetEdgeLength.hh>
#include <SurfaceMaps/Misc/ConstantCurvatureGeometry.hh>

namespace SurfaceMaps
{

/// This prevents triangles from degenerating or inverting or spanning a full hemisphere, and edges from approaching length pi
template <typename T>
T barrier_energy(
        const Vec3<T>& _a_sphere,
        const Vec3<T>& _b_sphere,
        const Vec3<T>& _c_sphere)
{
    // Inexact check
    const T volume = TinyAD::col_mat(_a_sphere, _b_sphere, _c_sphere).determinant();
    if (volume <= 0.0)
        return INFINITY;

    // Exact check
    if (!ccw_exclusive_3d(TinyAD::to_passive(_a_sphere), TinyAD::to_passive(_b_sphere), TinyAD::to_passive(_c_sphere)))
        return INFINITY;

    return -log(volume);
}

template double barrier_energy(const Vec3<double>&, const Vec3<double>&, const Vec3<double>&);
template TinyAD::Double<12,false> barrier_energy(const Vec3<TinyAD::Double<12,false>>&, const Vec3<TinyAD::Double<12,false>>&, const Vec3<TinyAD::Double<12,false>>&);
template TinyAD::Double<12,true> barrier_energy(const Vec3<TinyAD::Double<12,true>>&, const Vec3<TinyAD::Double<12,true>>&, const Vec3<TinyAD::Double<12,true>>&);
template TinyAD::Double<6,false> barrier_energy(const Vec3<TinyAD::Double<6,false>>&, const Vec3<TinyAD::Double<6,false>>&, const Vec3<TinyAD::Double<6,false>>&);
template TinyAD::Double<6,true> barrier_energy(const Vec3<TinyAD::Double<6,true>>&, const Vec3<TinyAD::Double<6,true>>&, const Vec3<TinyAD::Double<6,true>>&);

template <typename T>
T regularization_energy(
        const Vec3<T>& _a_sphere,
        const Vec3<T>& _b_sphere,
        const Vec3<T>& _c_sphere)
{
    T dirichlet = 0.0;

    dirichlet += (_a_sphere - _b_sphere).squaredNorm();
    dirichlet += (_b_sphere - _c_sphere).squaredNorm();
    dirichlet += (_c_sphere - _a_sphere).squaredNorm();

    return dirichlet;
}

template <typename T>
T map_energy(
        const Vec3<T>& _a_lifted_A,
        const Vec3<T>& _b_lifted_A,
        const Vec3<T>& _c_lifted_A,
        const Vec3<T>& _a_lifted_B,
        const Vec3<T>& _b_lifted_B,
        const Vec3<T>& _c_lifted_B)
{
    // Compute local 2D coordinate systems for lifted T triangles
    Vec2<T> a_lifted_local_A, b_lifted_local_A, c_lifted_local_A, a_lifted_local_B, b_lifted_local_B, c_lifted_local_B;
    to_local_coordinates(_a_lifted_A, _b_lifted_A, _c_lifted_A, a_lifted_local_A, b_lifted_local_A, c_lifted_local_A);
    to_local_coordinates(_a_lifted_B, _b_lifted_B, _c_lifted_B, a_lifted_local_B, b_lifted_local_B, c_lifted_local_B);

    // Matrices with edge vectors as columns
    Eigen::Matrix2<T> M_A;
    Eigen::Matrix2<T> M_B;
    M_A << b_lifted_local_A - a_lifted_local_A, c_lifted_local_A - a_lifted_local_A;
    M_B << b_lifted_local_B - a_lifted_local_B, c_lifted_local_B - a_lifted_local_B;

    // Compute areas of lifted triangles
    const T area_lifted_A = 0.5 * M_A.determinant();
    const T area_lifted_B = 0.5 * M_B.determinant();

    // Don't allow degenerate triangles
    if (area_lifted_A <= 0 || area_lifted_B <= 0)
        return INFINITY;

    // Compute map jacobian and its inverse
    Eigen::Matrix2<T> J = M_B * M_A.inverse();
    Eigen::Matrix2<T> J_inv = M_A * M_B.inverse();

    // Compute area-weighted symmetric Dirichlet energy
    return area_lifted_B * J.squaredNorm() + area_lifted_A * J_inv.squaredNorm();
}

/// Mesh energy based on symmetric Dirichlet energy to perfect equilateral triangles
template <typename T>
T mesh_energy(
        const Vec3<T>& _a_lifted_A,
        const Vec3<T>& _b_lifted_A,
        const Vec3<T>& _c_lifted_A,
        const Vec3<T>& _a_lifted_B,
        const Vec3<T>& _b_lifted_B,
        const Vec3<T>& _c_lifted_B,
        const T& _target_edge_length_A,
        const T& _target_edge_length_B)
{
    // Choose minimum of sizing picked on A and B
    const T tel = fmin(_target_edge_length_A, _target_edge_length_B);

    const Vec3<T> a_ideal = { 0.0, 0.0, 0.0 };
    const Vec3<T> b_ideal_A = { tel, 0.0, 0.0 };
    const Vec3<T> c_ideal_A = { 0.5 * tel, sqrt(3) * 0.5 * tel, 0 };
    const Vec3<T> b_ideal_B = { tel, 0.0, 0.0 };
    const Vec3<T> c_ideal_B = { 0.5 * tel, sqrt(3) * 0.5 * tel, 0 };

    T mesh_energy = 0.0;

    mesh_energy += map_energy(_a_lifted_A, _b_lifted_A, _c_lifted_A, a_ideal, b_ideal_A, c_ideal_A);
    mesh_energy += map_energy(_a_lifted_B, _b_lifted_B, _c_lifted_B, a_ideal, b_ideal_B, c_ideal_B);

    return mesh_energy;
}

/// Compute surface approx energy for all A or B vertices in a T triangle
template <typename T>
T surface_approx_energy(
        const Vec3<T>& _a_sphere,
        const Vec3<T>& _b_sphere,
        const Vec3<T>& _c_sphere,
        const Vec3<T>& _a_lifted,
        const Vec3<T>& _b_lifted,
        const Vec3<T>& _c_lifted,
        const std::vector<SVH>& _vhs_X,
        const TriMesh& _mesh_X,
        const ExternalProperty<VH, Vec3d>& _embedding_X,
        const ExternalProperty<VH, double>& _vertex_areas_X,
        const double _target_approx)
{
    T E_approx = 0.0;

    for (auto vh : _vhs_X)
    {
        // Compute barycentric coordinates (in sphere ambient space) of A/B vertex in T triangle
        T alpha, beta, gamma;
        barys_abc_3d(_a_sphere, _b_sphere, _c_sphere, _embedding_X[vh], Spherical, alpha, beta);
        gamma = 1.0 - alpha - beta;

        // Compute base point of A/B vertex on lifted T triangle
        Vec3<T> p_base = _a_lifted * alpha + _b_lifted * beta + _c_lifted * gamma;

        // Measure squared Euclidean distance of A/B vertex to base point on T triangle
        E_approx += (1.0 / sqr(_target_approx) * _vertex_areas_X[vh]) * (_mesh_X.point(vh) - p_base).squaredNorm();
    }

    return E_approx;
}

template <typename T>
T schueller_barrier(
        const T& _x, const double _s)
{
    if (_x <= 0.0)
    {
        return INFINITY;
    }
    else if (_x < _s)
    {
        const T g = _x * sqr(_x) * (1.0 / (_s * sqr(_s)))
            - 3.0 * sqr(_x) * (1.0 / sqr(_s))
            + 3.0 * _x * (1.0 / _s);
        return 1.0 / g - 1.0;
    }
    else
    {
        return 0.0;
    }
}

/// Compute surface approx energy for all A or B vertices in a T triangle
template <typename T>
T surface_approx_barrier(
        const Vec3<T>& _a_sphere,
        const Vec3<T>& _b_sphere,
        const Vec3<T>& _c_sphere,
        const Vec3<T>& _a_lifted,
        const Vec3<T>& _b_lifted,
        const Vec3<T>& _c_lifted,
        const std::vector<SVH>& _vhs_X,
        const TriMesh& _mesh_X,
        const ExternalProperty<VH, Vec3d>& _embedding_X,
        const ExternalProperty<VH, double>& _vertex_areas_X,
        const double& _bound)
{
    T E_approx_barrier = 0.0;

    for (auto vh : _vhs_X)
    {
        // Compute barycentric coordinates (in sphere ambient space) of A/B vertex in T triangle
        T alpha, beta, gamma;
        barys_abc_3d(_a_sphere, _b_sphere, _c_sphere, _embedding_X[vh], Spherical, alpha, beta);
        gamma = 1.0 - alpha - beta;

        // Compute base point of A/B vertex on lifted T triangle
        Vec3<T> p_base = _a_lifted * alpha + _b_lifted * beta + _c_lifted * gamma;

        // TODO: Use simplified formula from paper
        const T distance = (_mesh_X.point(vh) - p_base).norm();
        E_approx_barrier += _vertex_areas_X[vh] * schueller_barrier(_bound - distance, _bound);
    }

    return E_approx_barrier;
}

template double surface_approx_barrier(const Vec3<double>& _a_sphere, const Vec3<double>& _b_sphere, const Vec3<double>& _c_sphere,
        const Vec3<double>& _a_lifted, const Vec3<double>& _b_lifted, const Vec3<double>& _c_lifted,
        const std::vector<SVH>& _vhs_X, const TriMesh& _mesh_X, const ExternalProperty<VH, Vec3d>& _embedding_X,
        const ExternalProperty<VH, double>& _vertex_areas_X, const double& _bound);

template <typename T>
T angle_energy(
        const Vec3<T>& _a_lifted,
        const Vec3<T>& _b_lifted,
        const Vec3<T>& _c_lifted)
{
    const Vec3<T> ab = (_b_lifted - _a_lifted).normalized();
    const Vec3<T> bc = (_c_lifted - _b_lifted).normalized();
    const Vec3<T> ca = (_a_lifted - _c_lifted).normalized();

    const T alpha = acos(ab.dot(-ca));
    const T beta = acos(bc.dot(-ab));
    const T gamma = acos(ca.dot(-bc));

    return sqr(alpha - M_PI / 3.0)
         + sqr(beta - M_PI / 3.0)
         + sqr(gamma - M_PI / 3.0);
}

/// soft landmark energy
template <typename T>
T soft_landmark_energy(
        const Vec3d& _landmark_X,
        const Vec3<T>& _landmark_T)
{
    T dist_points = (_landmark_X - _landmark_T).squaredNorm();
    return dist_points;
}


/// Evaluate energy that only concern a single mesh for a triangle of mesh T (i.e., barrier, regulaization, surface approx)
template <typename T>
T eval_singlemesh_energy_triangle_T(
        const Vec3<T>& _a_sphere,
        const Vec3<T>& _b_sphere,
        const Vec3<T>& _c_sphere,
        const int& _mesh_idx,
        const MapState& _map_state,
        const std::vector<SVH>& _vhs_A, // A vertices in triangle of T
        const AdaptiveTriangulationsSettings& _settings)
{
    T E_barrier = 0.0;
    if (_settings.w_barrier > 0.0)
    {
        E_barrier += barrier_energy(_a_sphere, _b_sphere, _c_sphere);

        ISM_ASSERT_NOT_NAN(E_barrier);
        if (E_barrier == INFINITY)
            return INFINITY;
    }

    T E_reg = 0.0;
    if (_settings.w_reg > 0.0)
    {
        E_reg += regularization_energy(_a_sphere, _b_sphere, _c_sphere);

        ISM_ASSERT_NOT_NAN(E_reg);
        if (E_reg == INFINITY)
            return INFINITY;
    }

    T E_approx = 0.0;
    if (_settings.w_approx > 0.0)
    {
        const Vec3<T> a_lifted = lift_vertex_to_surface(_a_sphere, _map_state.meshes_input[_mesh_idx], _map_state.meshes_embeddings_input[_mesh_idx], _map_state.bsp_embeddings_input[_mesh_idx]);
        const Vec3<T> b_lifted = lift_vertex_to_surface(_b_sphere, _map_state.meshes_input[_mesh_idx], _map_state.meshes_embeddings_input[_mesh_idx], _map_state.bsp_embeddings_input[_mesh_idx]);
        const Vec3<T> c_lifted = lift_vertex_to_surface(_c_sphere, _map_state.meshes_input[_mesh_idx], _map_state.meshes_embeddings_input[_mesh_idx], _map_state.bsp_embeddings_input[_mesh_idx]);
        E_approx += surface_approx_energy(_a_sphere, _b_sphere, _c_sphere, a_lifted, b_lifted, c_lifted, _vhs_A, _map_state.meshes_input[_mesh_idx], _map_state.embeddings_input[_mesh_idx], _map_state.vertex_areas_input[_mesh_idx], _settings.approx_error);

        ISM_ASSERT_NOT_NAN(E_approx);
        if (E_approx == INFINITY)
            return INFINITY;
    }

    T E_approx_barrier = 0.0;
    if (_settings.w_bound > 0.0)
    {
        const Vec3<T> a_lifted = lift_vertex_to_surface(_a_sphere, _map_state.meshes_input[_mesh_idx], _map_state.meshes_embeddings_input[_mesh_idx], _map_state.bsp_embeddings_input[_mesh_idx]);
        const Vec3<T> b_lifted = lift_vertex_to_surface(_b_sphere, _map_state.meshes_input[_mesh_idx], _map_state.meshes_embeddings_input[_mesh_idx], _map_state.bsp_embeddings_input[_mesh_idx]);
        const Vec3<T> c_lifted = lift_vertex_to_surface(_c_sphere, _map_state.meshes_input[_mesh_idx], _map_state.meshes_embeddings_input[_mesh_idx], _map_state.bsp_embeddings_input[_mesh_idx]);
        E_approx_barrier += surface_approx_barrier(_a_sphere, _b_sphere, _c_sphere, a_lifted, b_lifted, c_lifted, _vhs_A, _map_state.meshes_input[_mesh_idx], _map_state.embeddings_input[_mesh_idx], _map_state.vertex_areas_input[_mesh_idx], _settings.surface_approx_bound);

        ISM_ASSERT_NOT_NAN(E_approx_barrier);
        if (E_approx_barrier == INFINITY)
            return INFINITY;
    }

    T E_angle = 0.0;
    if (_settings.w_angle > 0.0)
    {
        const Vec3<T> a_lifted = lift_vertex_to_surface(_a_sphere, _map_state.meshes_input[_mesh_idx], _map_state.meshes_embeddings_input[_mesh_idx], _map_state.bsp_embeddings_input[_mesh_idx]);
        const Vec3<T> b_lifted = lift_vertex_to_surface(_b_sphere, _map_state.meshes_input[_mesh_idx], _map_state.meshes_embeddings_input[_mesh_idx], _map_state.bsp_embeddings_input[_mesh_idx]);
        const Vec3<T> c_lifted = lift_vertex_to_surface(_c_sphere, _map_state.meshes_input[_mesh_idx], _map_state.meshes_embeddings_input[_mesh_idx], _map_state.bsp_embeddings_input[_mesh_idx]);
        E_angle += angle_energy(a_lifted, b_lifted, c_lifted);

        ISM_ASSERT_NOT_NAN(E_angle);
        if (E_angle == INFINITY)
            return INFINITY;
    }

    return (E_barrier * _settings.w_barrier
            + E_reg * _settings.w_reg
            + E_approx * _settings.w_approx
            + E_approx_barrier * _settings.w_bound
            + E_angle * _settings.w_angle) * (1.0/(double)_map_state.meshes_input.size());

}

// Explicit template instantiation
template double eval_singlemesh_energy_triangle_T(const Vec3<double>&, const Vec3<double>&, const Vec3<double>&, const int&, const MapState&, const std::vector<SVH>&, const AdaptiveTriangulationsSettings&);
template TinyAD::Double<6, false> eval_singlemesh_energy_triangle_T(const Vec3<TinyAD::Double<6, false>>&, const Vec3<TinyAD::Double<6, false>>&, const Vec3<TinyAD::Double<6, false>>&, const int&, const MapState&, const std::vector<SVH>&, const AdaptiveTriangulationsSettings&);
template TinyAD::Double<6, true> eval_singlemesh_energy_triangle_T(const Vec3<TinyAD::Double<6, true>>&, const Vec3<TinyAD::Double<6, true>>&, const Vec3<TinyAD::Double<6, true>>&, const int&, const MapState&, const std::vector<SVH>&, const AdaptiveTriangulationsSettings&);


/// Evaluate energy that only pair of meshes for a triangle of mesh T (i.e., map energy)
template <typename T>
T eval_trianglepair_energy_triangle_T(
        const Vec3<T>& _a_sphere_A,
        const Vec3<T>& _b_sphere_A,
        const Vec3<T>& _c_sphere_A,
        const Vec3<T>& _a_sphere_B,
        const Vec3<T>& _b_sphere_B,
        const Vec3<T>& _c_sphere_B,
        const int& _mesh_A_idx,
        const int& _mesh_B_idx,
        const MapState& _map_state,
        const AdaptiveTriangulationsSettings& _settings)
{
    if (_map_state.pairs_map_distortion.size() == 0)
        return 0.0;

    // Compute lifted positions of vertices
    const Vec3<T> a_lifted_A = lift_vertex_to_surface(_a_sphere_A, _map_state.meshes_input[_mesh_A_idx], _map_state.meshes_embeddings_input[_mesh_A_idx], _map_state.bsp_embeddings_input[_mesh_A_idx]);
    const Vec3<T> b_lifted_A = lift_vertex_to_surface(_b_sphere_A, _map_state.meshes_input[_mesh_A_idx], _map_state.meshes_embeddings_input[_mesh_A_idx], _map_state.bsp_embeddings_input[_mesh_A_idx]);
    const Vec3<T> c_lifted_A = lift_vertex_to_surface(_c_sphere_A, _map_state.meshes_input[_mesh_A_idx], _map_state.meshes_embeddings_input[_mesh_A_idx], _map_state.bsp_embeddings_input[_mesh_A_idx]);
    const Vec3<T> a_lifted_B = lift_vertex_to_surface(_a_sphere_B, _map_state.meshes_input[_mesh_B_idx], _map_state.meshes_embeddings_input[_mesh_B_idx], _map_state.bsp_embeddings_input[_mesh_B_idx]);
    const Vec3<T> b_lifted_B = lift_vertex_to_surface(_b_sphere_B, _map_state.meshes_input[_mesh_B_idx], _map_state.meshes_embeddings_input[_mesh_B_idx], _map_state.bsp_embeddings_input[_mesh_B_idx]);
    const Vec3<T> c_lifted_B = lift_vertex_to_surface(_c_sphere_B, _map_state.meshes_input[_mesh_B_idx], _map_state.meshes_embeddings_input[_mesh_B_idx], _map_state.bsp_embeddings_input[_mesh_B_idx]);

    T E_map = 0.0;
    if (_settings.w_map > 0.0)
    {
        E_map += map_energy(a_lifted_A, b_lifted_A, c_lifted_A, a_lifted_B, b_lifted_B, c_lifted_B);

        ISM_ASSERT_NOT_NAN(E_map);
        if (E_map == INFINITY)
            return INFINITY;
    }
    return E_map * _settings.w_map * (1.0/(double)_map_state.pairs_map_distortion.size());  // Weight map energy like other energies
}
// Explicit template instantiation
template double eval_trianglepair_energy_triangle_T(const Vec3<double>&, const Vec3<double>&, const Vec3<double>&, const Vec3<double>&, const Vec3<double>&, const Vec3<double>&, const int&, const int&, const MapState&, const AdaptiveTriangulationsSettings&);
template TinyAD::Double<12, false> eval_trianglepair_energy_triangle_T(const Vec3<TinyAD::Double<12, false>>&, const Vec3<TinyAD::Double<12, false>>&, const Vec3<TinyAD::Double<12, false>>&, const Vec3<TinyAD::Double<12, false>>&, const Vec3<TinyAD::Double<12, false>>&, const Vec3<TinyAD::Double<12, false>>&, const int&, const int&, const MapState&, const AdaptiveTriangulationsSettings&);
template TinyAD::Double<12, true> eval_trianglepair_energy_triangle_T(const Vec3<TinyAD::Double<12, true>>&, const Vec3<TinyAD::Double<12, true>>&, const Vec3<TinyAD::Double<12, true>>&, const Vec3<TinyAD::Double<12, true>>&, const Vec3<TinyAD::Double<12, true>>&, const Vec3<TinyAD::Double<12, true>>&, const int&, const int&, const MapState&, const AdaptiveTriangulationsSettings&);


/// Evaluate mesh energy for triangle of T for MapState
template <typename T>
T eval_mesh_energy_triangle_T(
        const Vec3<T>& _a_sphere,
        const Vec3<T>& _b_sphere,
        const Vec3<T>& _c_sphere,
        const int& _mesh_idx,
        const T& _tel,
        const MapState& _map_state,
        const AdaptiveTriangulationsSettings& _settings)
{
    T E_mesh = 0.0;
    if (_settings.w_mesh > 0.0)
    {
        const Vec3<T> a_lifted = lift_vertex_to_surface(_a_sphere, _map_state.meshes_input[_mesh_idx], _map_state.meshes_embeddings_input[_mesh_idx], _map_state.bsp_embeddings_input[_mesh_idx]);
        const Vec3<T> b_lifted = lift_vertex_to_surface(_b_sphere, _map_state.meshes_input[_mesh_idx], _map_state.meshes_embeddings_input[_mesh_idx], _map_state.bsp_embeddings_input[_mesh_idx]);
        const Vec3<T> c_lifted = lift_vertex_to_surface(_c_sphere, _map_state.meshes_input[_mesh_idx], _map_state.meshes_embeddings_input[_mesh_idx], _map_state.bsp_embeddings_input[_mesh_idx]);

        const Vec3<T> a_ideal = { 0.0, 0.0, 0.0 };
        const Vec3<T> b_ideal_A = { _tel, 0.0, 0.0 };
        const Vec3<T> c_ideal_A = { 0.5 * _tel, sqrt(3) * 0.5 * _tel, 0 };

        E_mesh += map_energy(a_lifted, b_lifted, c_lifted, a_ideal, b_ideal_A, c_ideal_A);

        ISM_ASSERT_NOT_NAN(E_mesh);
        if (E_mesh == INFINITY)
            return INFINITY;
    }

    return (E_mesh * _settings.w_mesh) * (1.0/(double)_map_state.meshes_input.size());
}
// Explicit template instantiation
template double eval_mesh_energy_triangle_T(const Vec3<double>&, const Vec3<double>&, const Vec3<double>&, const int&, const double&, const MapState&, const AdaptiveTriangulationsSettings&);
template TinyAD::Double<12, false> eval_mesh_energy_triangle_T(const Vec3<TinyAD::Double<12, false>>&, const Vec3<TinyAD::Double<12, false>>&, const Vec3<TinyAD::Double<12, false>>&, const int&, const TinyAD::Double<12, false>&, const MapState&, const AdaptiveTriangulationsSettings&);
template TinyAD::Double<12, true> eval_mesh_energy_triangle_T(const Vec3<TinyAD::Double<12, true>>&, const Vec3<TinyAD::Double<12, true>>&, const Vec3<TinyAD::Double<12, true>>&, const int&, const TinyAD::Double<12, true>&, const MapState&, const AdaptiveTriangulationsSettings&);

/// Evaluate energy for a landmark vertex of mesh T
template <typename T>
T eval_landmark_T_energy(
        const Vec3<T>& _p_sphere,
        const VH& _vh_input,
        const int& _idx_input,
        const MapState& _map_state,
        const AdaptiveTriangulationsSettings& _settings)
{
    T E_landmark_world = 0.0;
    if (_settings.w_landmark_world > 0.0)
    {
        Vec3<T> p_lifted_A = lift_vertex_to_surface(_p_sphere, _map_state.meshes_input[_idx_input], _map_state.meshes_embeddings_input[_idx_input], _map_state.bsp_embeddings_input[_idx_input]);
        E_landmark_world += soft_landmark_energy(_map_state.meshes_input[_idx_input].point(_vh_input), p_lifted_A);

        ISM_ASSERT_NOT_NAN(E_landmark_world);
    }

    T E_landmark_sphere = 0.0;
    if (_settings.w_landmark_sphere > 0.0)
    {
        E_landmark_sphere += soft_landmark_energy(_map_state.embeddings_input[_idx_input][_vh_input], _p_sphere);

        ISM_ASSERT_NOT_NAN(E_landmark_sphere);
    }

    return (_settings.w_landmark_world * E_landmark_world
         + _settings.w_landmark_sphere * E_landmark_sphere) * (1.0/(double)_map_state.meshes_input.size() /(double)_map_state.landmarks_T.size());
}

// Explicit template instantiation
template double eval_landmark_T_energy(const Vec3<double>&, const VH&,  const int&, const MapState&, const AdaptiveTriangulationsSettings&);
template TinyAD::Double<2, false> eval_landmark_T_energy(const Vec3<TinyAD::Double<2, false>>&, const VH&,  const int&, const MapState&, const AdaptiveTriangulationsSettings&);
template TinyAD::Double<2, true> eval_landmark_T_energy(const Vec3<TinyAD::Double<2, true>>&, const VH&,  const int&, const MapState&, const AdaptiveTriangulationsSettings&);

}
