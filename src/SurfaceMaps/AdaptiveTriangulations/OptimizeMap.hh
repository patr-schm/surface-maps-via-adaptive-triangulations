/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Dörte Pieper, Patrick Schmidt
 */
#pragma once

#include <SurfaceMaps/AdaptiveTriangulations/MapState.hh>
#include <SurfaceMaps/AdaptiveTriangulations/AdaptiveTriangulationsSettings.hh>

#include <SurfaceMaps/Utils/Filesystem.hh>

namespace SurfaceMaps
{

/**
 * To return for outer loops for stopping criteria
 */
enum OptimizeMapStatus
{
    normal,
    under_threshold,
    line_search_converge
};

/**
 * Continuous manifold optimization of _map_state
 */
OptimizeMapStatus optimize_map(
        MapState& _map_state,
        const int& _iterations,
        std::function<void(const std::string&)> _callback = [] (const std::string&) {},
        const fs::path& _csv_path = "",
        const AdaptiveTriangulationsSettings& _settings = AdaptiveTriangulationsSettings());

}
