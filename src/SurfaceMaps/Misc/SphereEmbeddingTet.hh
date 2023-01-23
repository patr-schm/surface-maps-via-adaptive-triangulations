/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt
 */
#pragma once

#include <SurfaceMaps/Types.hh>

namespace SurfaceMaps
{

/// Simple sphere embedding of genus 0 mesh.
/// (1) Pick two neighboring triangles (a, b, c), (c, a, d) and embed them as two faces of a regular tetrahedron.
/// (2) Find a path between the opposite vertices b, d and embed it via constant speed parametrization.
/// (3) Embed all other vertices on the surface of the tetrahedron via harmonic parametrization and project them to the sphere.
bool sphere_embedding_tet(
        const TriMesh& _mesh_L,
        const EH _eh_main,
        ExternalProperty<VH, Vec3d>& _embedding);

}
