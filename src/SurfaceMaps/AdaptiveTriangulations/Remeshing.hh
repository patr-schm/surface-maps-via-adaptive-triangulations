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
#include <SurfaceMaps/AdaptiveTriangulations/AdaptiveTriangulationsSettings.hh>

#include <SurfaceMaps/Utils/Filesystem.hh>

namespace SurfaceMaps
{

/**
 * One iteration of remeshing of mesh_T and its embeddings T_A and T_B.
 * Return true if some operations change the map_state.
 */
bool remesh_T(
        MapState& _map_state,
        const AdaptiveTriangulationsSettings& _settings,
        const fs::path& _csv_path = "",
        std::function<void (const std::string&)> _callback = [] (const std::string&) {});

}
