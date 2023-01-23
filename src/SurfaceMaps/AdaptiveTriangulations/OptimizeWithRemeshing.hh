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
 * Maximum distance of a landmark to its target position in world space.
 */
double landmark_error(
        const MapState& _map_state);

/**
 * Initialize target edge lengths fields for mesh_A and mesh_B in MapState.
 */
void init_target_edge_lengths(
        MapState& _map_state,
        const AdaptiveTriangulationsSettings& _settings);

/**
  * Takes sizing fields as input.
  * Constructs bsp trees.
  */
void optimize_with_remeshing(
        MapState& _map_state,
        const std::vector<ExternalProperty<VH, double>>& _mesh_tels,
        const AdaptiveTriangulationsSettings& _settings = AdaptiveTriangulationsSettings(),
        const fs::path& _csv_path = "",
        std::function<void (const std::string&)> _callback_for_optim = [] (const std::string&) {},
        std::function<void (const std::string&)> _callback_for_remesh = [] (const std::string&) {});

/**
  * Computes sizing fields before calling main algorithm.
  * Constructs bsp trees.
  */
void optimize_with_remeshing(
        MapState& _map_state,
        const AdaptiveTriangulationsSettings& _settings = AdaptiveTriangulationsSettings(),
        const fs::path& _csv_path = "",
        std::function<void (const std::string&)> _callback_for_optim = [] (const std::string&) {},
        std::function<void (const std::string&)> _callback_for_remesh = [] (const std::string&) {});

}
