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
#include <SurfaceMaps/AdaptiveTriangulations/MapState.hh>

namespace SurfaceMaps
{

// Return TriMesh of sphere embedding
TriMesh embedding_to_mesh(
        const TriMesh& _mesh_orig,
        const ExternalProperty<VH, Vec3d>& _embedding);

std::vector<TriMesh> lifted_meshes_from_mapstate(
        const MapState& _map_state);

/// Return an orthogonal vector to input vector
Vec3d any_orthogonal(
        const Vec3d& _v);

void compute_local_bases(
        const ExternalProperty<VH, Vec3d>& _embedding_input,
        ExternalProperty<VH, Vec3d>& _B1, // maxbe external property?
        ExternalProperty<VH, Vec3d>& _B2,
        const TriMesh& _mesh_input);

/// Compute local coordinates for input points, coordinate system centered in a
template <typename T>
void to_local_coordinates(
        const Vec3<T>& _a_in,
        const Vec3<T>& _b_in,
        const Vec3<T>& _c_in,
        Vec2<T>& _a_out,
        Vec2<T>& _b_out,
        Vec2<T>& _c_out)
{
    Vec3<T> normal = ((_b_in - _a_in).cross(_c_in - _a_in)).normalized();
    Vec3<T> basis0 = (_b_in - _a_in).normalized();
    Vec3<T> basis1 = normal.cross(basis0);

    _a_out = Vec2<T>(0.0, 0.0);
    _b_out = Vec2<T>((_b_in - _a_in).dot(basis0), 0.0);
    _c_out = Vec2<T>((_c_in - _a_in).dot(basis0), (_c_in - _a_in).dot(basis1));
}

/// Return number of flipped triangles on sphere with volume calculation
int n_flipped_volume(const MapState& _map_state);

/// compute length of the diagonal of the axis aligned bounding box of the mesh
double length_diagonal_of_bounding_box(const TriMesh& _mesh);

}
