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

namespace SurfaceMaps
{

/// Compute colors per face for distortion between zwo meshes
ExternalProperty<FH, Color> compute_distortion_heatmap(
        const TriMesh& _mesh_1,
        const TriMesh& _mesh_2,
        const float &_range_max = 100.0f);

}

