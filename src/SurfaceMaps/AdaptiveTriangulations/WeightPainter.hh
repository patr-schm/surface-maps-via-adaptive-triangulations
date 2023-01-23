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
#include <SurfaceMaps/Utils/Filesystem.hh>

namespace SurfaceMaps
{

/// Viewer with UI for painting and saving a scalar field for a mesh
void view_weight_painter(
        const TriMesh& _mesh,
        ExternalProperty<VH, double>& _tels,
        const fs::path& _target_dir_scalar = "");

}
