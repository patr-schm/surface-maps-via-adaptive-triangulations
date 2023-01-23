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

namespace SurfaceMaps
{

constexpr double default_coarse_approx_error = 0.01;
constexpr double default_fine_approx_error = 0.001;
constexpr double default_w_approx = 1.0;

AdaptiveTriangulationsSettings landmark_phase_settings();

AdaptiveTriangulationsSettings coarse_phase_settings(
        const double _target_approximation_error = default_coarse_approx_error,
        const double _w_approx = default_w_approx);

AdaptiveTriangulationsSettings fine_phase_settings(
        const double _target_approximation_error = default_fine_approx_error,
        const double _w_approx = default_w_approx);

void landmark_phase(
        MapState& _map_state,
        std::function<void (const std::string&)> _iteration_callback = [] (const std::string&) {});

void coarse_phase(
        MapState& _map_state,
        const double _target_approximation_error = default_coarse_approx_error,
        const double _w_approx = default_w_approx,
        std::function<void (const std::string&)> _iteration_callback = [] (const std::string&) {});

void fine_phase(
        MapState& _map_state,
        const double _target_approximation_error = default_fine_approx_error,
        const double _w_approx = default_w_approx,
        std::function<void (const std::string&)> _iteration_callback = [] (const std::string&) {});

void release_landmarks(
        MapState& _map_state);

void satisfy_approx_bound_local_refinement(
        MapState& _map_state,
        const double &_bound,
        const int& _max_iters = 20,
        std::function<void (const std::string&)> _iteration_callback = [] (auto) {});

void optimize_coarse_to_fine(
        MapState& _map_state,
        std::function<void (const std::string&)> _iteration_callback = [] (const std::string&) {},
        std::function<void (const std::string&)> _phase_callback = [] (const std::string&) {});

}
