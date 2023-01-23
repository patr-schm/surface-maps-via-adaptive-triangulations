/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Dörte Pieper, Patrick Schmidt
 */
#include <SurfaceMaps/AdaptiveTriangulations/OptimizeCoarseToFine.hh>

#include <SurfaceMaps/AdaptiveTriangulations/OptimizeWithRemeshing.hh>
#include <SurfaceMaps/AdaptiveTriangulations/GuaranteeSurfaceApproximation.hh>
#include <SurfaceMaps/AdaptiveTriangulations/AssignVerticesToFaces.hh>
#include <TinyAD/Utils/Timer.hh>

namespace SurfaceMaps
{

AdaptiveTriangulationsSettings landmark_phase_settings()
{
    AdaptiveTriangulationsSettings settings;

    settings.max_iterations = 100;
    settings.set_adaptive_tel_params(0.01, 0.001, 100.0);
    settings.w_barrier = 1.0;
    settings.w_reg = 0.0;
    settings.w_map = 0.0;
    settings.w_mesh = 0.0;
    settings.w_approx = 0.0;
    settings.hard_landmark_constraints = false;
    settings.w_landmark_world = 0.0;
    settings.w_landmark_sphere = 1e6;
    settings.conv_thresh_newton_decrement = 0.0;
    settings.landmark_threshold = 1e-4;

    return settings;
}

AdaptiveTriangulationsSettings coarse_phase_settings(
        const double _target_approximation_error,
        const double _w_approx)
{
    AdaptiveTriangulationsSettings settings;

    settings.max_iterations = 50;
    settings.set_adaptive_tel_params(_target_approximation_error, 0.001, 100.0);
    settings.w_barrier = 1e-3;
    settings.w_reg = 0.0;
    settings.w_map = 1.0;
    settings.w_mesh = 1.0;
    settings.w_approx = _w_approx;
    settings.hard_landmark_constraints = true;
    settings.w_landmark_world = 0.0;
    settings.w_landmark_sphere = 0.0;
    settings.conv_thresh_newton_decrement = 1e-4;
    settings.landmark_threshold = -INF_DOUBLE;

    return settings;
}

AdaptiveTriangulationsSettings fine_phase_settings(
        const double _target_approximation_error,
        const double _w_approx)
{
    AdaptiveTriangulationsSettings settings;

    settings.max_iterations = 50;
    settings.set_adaptive_tel_params(_target_approximation_error, 0.001, 100.0);
    settings.w_barrier = 1e-6;
    settings.w_reg = 0.0;
    settings.w_map = 1.0;
    settings.w_mesh = 1.0;
    settings.w_approx = _w_approx;
    settings.hard_landmark_constraints = true;
    settings.w_landmark_world = 0.0;
    settings.w_landmark_sphere = 0.0;
    settings.conv_thresh_newton_decrement = 1e-4;
    settings.landmark_threshold = -INF_DOUBLE;

    return settings;
}

void landmark_phase(
        MapState& _map_state,
        std::function<void (const std::string&)> _iteration_callback)
{
    ISM_INFO("#################### Starting Landmark Phase ####################");

    AdaptiveTriangulationsSettings settings = landmark_phase_settings();
    optimize_with_remeshing(_map_state, settings, "", _iteration_callback, _iteration_callback);

    ISM_DEBUG_VAR(landmark_error(_map_state));
    ISM_EXPECT_L(landmark_error(_map_state), settings.landmark_threshold);
}

void coarse_phase(
        MapState& _map_state,
        const double _target_approximation_error,
        const double _w_approx,
        std::function<void (const std::string&)> _iteration_callback)
{
    ISM_INFO("#################### Starting Coarse Phase ####################");

    AdaptiveTriangulationsSettings settings = coarse_phase_settings(_target_approximation_error, _w_approx);
    optimize_with_remeshing(_map_state, settings, "", _iteration_callback, _iteration_callback);
}

void fine_phase(
        MapState& _map_state,
        const double _target_approximation_error,
        const double _w_approx,
        std::function<void (const std::string&)> _iteration_callback)
{
    ISM_INFO("#################### Starting Fine Phase ####################");

    AdaptiveTriangulationsSettings settings = fine_phase_settings(_target_approximation_error, _w_approx);
    optimize_with_remeshing(_map_state, settings, "", _iteration_callback, _iteration_callback);
}

void release_landmarks(
        MapState& _map_state)
{
    _map_state.landmarks_T.clear();
    for(int i = 0; i < (int)_map_state.meshes_input.size(); ++i)
        _map_state.landmarks_input[i].clear();
}

void satisfy_approx_bound_local_refinement(
        MapState& _map_state,
        const double &_bound,
        const int& _max_iters,
        std::function<void (const std::string&)> _iteration_callback)
{
    // Re-compute vertex A/B to face T assignments
    assign_vertices_to_T_faces(_map_state);

    const double approx_error_factor = 0.33;
    AdaptiveTriangulationsSettings bound_satisfy_settings = fine_phase_settings(approx_error_factor * _bound);
    bound_satisfy_settings.max_iterations = 2;
    bound_satisfy_settings.optimize_per_remeshing_iters = 0;
    bound_satisfy_settings.w_map = 0.0;

    if (n_verts_over_bound(_map_state, _bound) == 0)
    {
        ISM_INFO("Surface approximation bound already satisfied.");
        return;
    }

    for (int i = 0; i < _max_iters; i++)
    {
        ISM_DEBUG_OUT("Satisfying surface approximation bound. Iteration: " << i);
        ISM_DEBUG_OUT("Current max distance: " << max_dist_verts_T(_map_state));

        ISM_DEBUG_VAR(max_dist_verts_T(_map_state));
        ISM_DEBUG_VAR(n_verts_over_bound(_map_state, _bound));

        local_refine_tel(_map_state, _bound, 0.5, true);
        optimize_with_remeshing(_map_state, _map_state.tels_input, bound_satisfy_settings, "", _iteration_callback, _iteration_callback);

        if (n_verts_over_bound(_map_state, _bound) == 0)
        {
            ISM_INFO("Bound satisfied in iteration " << i);
            ISM_DEBUG_OUT("Vertex num: " << _map_state.mesh_T.n_vertices());
            break;
        }

    }

    ISM_DEBUG_VAR(max_dist_verts_T(_map_state));
    ISM_ASSERT_EQ(n_verts_over_bound(_map_state, _bound), 0);
}

void optimize_coarse_to_fine(
        MapState& _map_state,
        std::function<void (const std::string&)> _iteration_callback,
        std::function<void (const std::string&)> _phase_callback)
{
    _phase_callback("init");

    TinyAD::Timer timer_landmark_phase("Landmark phase");
    landmark_phase(_map_state, _iteration_callback);
    timer_landmark_phase.stop();

    _phase_callback("after_landmarks");

    TinyAD::Timer timer_coarse_phase("Coarse phase");
    coarse_phase(_map_state, default_coarse_approx_error, default_w_approx, _iteration_callback);
    timer_coarse_phase.stop();

    _phase_callback("after_coarse");

    TinyAD::Timer timer_fine_phase("Fine phase");
    fine_phase(_map_state, default_fine_approx_error, default_w_approx, _iteration_callback);
    timer_fine_phase.stop();

    _phase_callback("after_fine");
}

}
