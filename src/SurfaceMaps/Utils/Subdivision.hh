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
#include <SurfaceMaps/Utils/RefinementMap.hh>

namespace SurfaceMaps
{

/// Subdivides mesh by inserting one new vertex per edge midpoint.
/// Returns a RefinementMap which allows to transfer properties (e.g. uv
/// coordinates) from the original to the refined mesh.
RefinementMap
subdivide_1_to_4(
        TriMesh& _mesh);

/// Subdivides mesh by inserting one new vertex per edge midpoint.
/// Expresses the resulting RefinementMap w.r.t. another mapping
/// _orig_ref_map from a previous refinement.
/// Use for iterated subdivision.
RefinementMap
subdivide_1_to_4(
        TriMesh& _mesh,
        const RefinementMap& _orig_ref_map);

}
