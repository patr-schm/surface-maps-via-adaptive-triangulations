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
#include <SurfaceMaps/AdaptiveTriangulations/LiftToSurface.hh>
#include <SurfaceMaps/Utils/BSPTree.hh>

namespace SurfaceMaps
{

/// compute adaptive target edge length for a mesh
ExternalProperty<VH, double> target_edge_lengths(
        const TriMesh& _mesh,
        const double _approx_error,
        const double _min_edge_length,
        const double _max_edge_length);

/// Pick target edge length for single vertex of T (for one _mesh_input)
template <typename T>
T pick_target_edge_length(
        const Vec3<T>& _p_sphere,
        const TriMesh& _mesh_input, // mesh_A or mesh_B
        const TriMesh& _mesh_input_embedding,
        const ExternalProperty<VH, double>& _mesh_input_tels,
        const BSPTree& _bsp_input)
{
    SFH f;
    T alpha, beta, gamma;
    bsp_tree_barys_face(_p_sphere, _mesh_input_embedding, _bsp_input, alpha, beta, gamma, f);

    // Interpolate p/w lin scalar field
    SVH vh_a, vh_b, vh_c;
    handles(_mesh_input, f, vh_a, vh_b, vh_c);
    T interpol_tel = _mesh_input_tels[vh_a] * alpha + _mesh_input_tels[vh_b] * beta + _mesh_input_tels[vh_c] * gamma;

    return interpol_tel;
}

/// Target edge length for single face from its three vertex target edge lengths
template <typename T>
T triangle_tel_from_vertex_tels(
        const T& _tel_v1,
        const T& _tel_v2,
        const T& _tel_v3,
        const bool _pick_minimum)
{
    if (_pick_minimum)
        return fmin(_tel_v1, fmin(_tel_v2, _tel_v3));
    else
        return (_tel_v1 + _tel_v2 + _tel_v3) / 3.0; // avg of vertices tel
}

/// Target edge length for triangle of T (on one _mesh_input)
template <typename T>
T target_edge_length(
        const Vec3<T>& _v0_sphere,
        const Vec3<T>& _v1_sphere,
        const TriMesh& _mesh_input, // mesh_A or mesh_B
        const TriMesh& _mesh_input_embedding,
        const ExternalProperty<VH, double>& _mesh_input_tels,
        const BSPTree& _bsp_input)
{
    Vec3<T> p = (_v0_sphere + _v1_sphere ) * 0.5; // Skip normalization because only used for barycentric coordinates
    return pick_target_edge_length(p, _mesh_input, _mesh_input_embedding, _mesh_input_tels, _bsp_input);
}

/// Target edge length for triangle of T (on one _mesh_input)
template <typename T>
T ideal_triangle_edge_length(
        const Vec3<T>& _a_sphere,
        const Vec3<T>& _b_sphere,
        const Vec3<T>& _c_sphere,
        const TriMesh& _mesh_input, // mesh_A or mesh_B
        const TriMesh& _mesh_input_embedding,
        const ExternalProperty<VH, double>& _mesh_input_tels,
        const BSPTree& _bsp_input)
{
    Vec3<T> middle_p = (_a_sphere + _b_sphere + _c_sphere) * (1.0/3.0); // Skip normalization because only used for barycentric coordinates
    return pick_target_edge_length(middle_p, _mesh_input, _mesh_input_embedding, _mesh_input_tels, _bsp_input);
}

}
