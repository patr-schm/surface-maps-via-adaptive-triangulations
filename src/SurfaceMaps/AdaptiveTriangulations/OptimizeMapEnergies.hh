/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Dörte Pieper, Patrick Schmidt
 */
#pragma once

#include <SurfaceMaps/Types.hh>
#include <SurfaceMaps/AdaptiveTriangulations/AdaptiveTriangulationsSettings.hh>
#include <SurfaceMaps/AdaptiveTriangulations/MapState.hh>

namespace SurfaceMaps
{

/// This prevents triangles from degenerating or inverting or spanning a full hemisphere, and edges from approaching length pi
template <typename T>
T barrier_energy(
        const Vec3<T>& _a_sphere,
        const Vec3<T>& _b_sphere,
        const Vec3<T>& _c_sphere);

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
        const double& _bound);

/// Evaluate energy that only concern a single mesh for a triangle of mesh T (i.e., barrier, regulaization, surface approx)
template <typename T>
T eval_singlemesh_energy_triangle_T(
        const Vec3<T>& _a_sphere,
        const Vec3<T>& _b_sphere,
        const Vec3<T>& _c_sphere,
        const int& _mesh_idx,
        const MapState& _map_state,
        const std::vector<SVH>& _vhs_A, // A vertices in triangle of T
        const AdaptiveTriangulationsSettings& _settings);

/// Evaluate energy that concerns pair of triangles on two meshes for a triangle of mesh T (i.e., map energy)
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
        const AdaptiveTriangulationsSettings& _settings);

/// Evaluate mesh energy for triangle of T for MapState
template <typename T>
T eval_mesh_energy_triangle_T(
        const Vec3<T>& _a_sphere,
        const Vec3<T>& _b_sphere,
        const Vec3<T>& _c_sphere,
        const int& _mesh_idx,
        const T& _tel,
        const MapState& _map_state,
        const AdaptiveTriangulationsSettings& _settings);

/// Evaluate energy for a landmark vertex of mesh T with two embeddings
template <typename T>
T eval_landmark_T_energy(
        const Vec3<T>& _p_sphere,
        const VH& _vh_input,
        const int& _idx_input,
        const MapState& _map_state,
        const AdaptiveTriangulationsSettings& _settings);

}
