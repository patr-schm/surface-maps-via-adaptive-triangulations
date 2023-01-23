/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Dörte Pieper, Patrick Schmidt
 */
#include <SurfaceMaps/AdaptiveTriangulations/OptimizeWithRemeshing.hh>

#include <SurfaceMaps/AdaptiveTriangulations/AdaptiveTriangulationsSettings.hh>
#include <SurfaceMaps/AdaptiveTriangulations/LiftToSurface.hh>
#include <SurfaceMaps/AdaptiveTriangulations/Helpers.hh>
#include <SurfaceMaps/AdaptiveTriangulations/Remeshing.hh>
#include <SurfaceMaps/AdaptiveTriangulations/OptimizeMap.hh>
#include <SurfaceMaps/AdaptiveTriangulations/AdaptiveTargetEdgeLength.hh>
#include <SurfaceMaps/AdaptiveTriangulations/AssignVerticesToFaces.hh>
#include <SurfaceMaps/AdaptiveTriangulations/OptimizeMapEnergies.hh>

#include <TinyAD/Utils/Timer.hh>

namespace SurfaceMaps
{

double landmark_error(
        const MapState& _map_state)
{
    for(int i = 0; i < (int)_map_state.meshes_input.size(); i++)
    {
        ISM_ASSERT_EQ(_map_state.landmarks_input[i].size(), _map_state.landmarks_T.size());
    }

    double max_error = 0.0;
    for (int i = 0; i < (int)_map_state.landmarks_T.size(); ++i)
    {
        const VH v_T = _map_state.landmarks_T[i];
        for(int j = 0; j < (int)_map_state.meshes_input.size(); j++)
        {
            const VH v_j = _map_state.landmarks_input[j][i];
            const Vec3d p_T_j = _map_state.embeddings_T[j][v_T];
            const Vec3d p_j = lift_vertex_to_surface(p_T_j, _map_state.meshes_input[j], _map_state.meshes_embeddings_input[j], _map_state.bsp_embeddings_input[j]);
            const Vec3d p_j_target = _map_state.meshes_input[j].point(v_j);
            // Euclidean distance between current and target position in world space
            max_error = std::max(max_error, (p_j - p_j_target).norm());
        }
    }

    return max_error;
}



void init_target_edge_lengths(
        MapState& _map_state,
        const AdaptiveTriangulationsSettings& _settings)
{
    _map_state.tels_input.clear();

    // Input for A, compute factor for other meshes based on mesh areas
    std::vector<double> factor_A_to_i_area;
    for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
        factor_A_to_i_area.push_back(total_area(_map_state.meshes_input[i]) / total_area(_map_state.meshes_input[0]));

    for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
        _map_state.tels_input.push_back(target_edge_lengths(_map_state.meshes_input[i], _settings.approx_error * sqrt(factor_A_to_i_area[i]), _settings.min_edge_length * sqrt(factor_A_to_i_area[i]), _settings.max_edge_length * sqrt(factor_A_to_i_area[i])));

    ISM_ASSERT_EQ(_map_state.meshes_input.size(), _map_state.tels_input.size());
}


namespace
{

void optimize_with_remeshing_impl(
        MapState& _map_state,
        const AdaptiveTriangulationsSettings& _settings,
        const fs::path& _csv_path,
        std::function<void (const std::string&)> _callback_for_optim,
        std::function<void (const std::string&)> _callback_for_remesh)
{
    const int n_meshes = _map_state.meshes_input.size();

//    // Init CSV file
//    if (_cur_iter == 0 && !_csv_path.empty())
//    {
//        make_file_directory(_csv_path);
//        fs::remove(_csv_path);
//        if (!_settings.energies_in_csv)
//            append_to_csv(_csv_path,
//                          "iteration", "iter_type", "derivatives (seconds)", "solve (seconds)", "line search (seconds)", "iteration (seconds)", "objective",
//                          "splits (seconds)", "n_splits", "collapses (seconds)", "n_collapses", "flips (seconds)", "n_flips", "n_vertices");
//        else
//            ISM_ERROR_throw("Fixme");
//    }

    for (int i = 0; i < n_meshes; i++)
        ISM_ASSERT(_map_state.tels_input[i].size_okay(_map_state.meshes_input[i]));

    // Assign to each T triangle the vertices of input meshes that are inside
    if (_settings.w_approx > 0.0)
    {
        assign_vertices_to_T_faces(_map_state);
    }
    else
    {
        _map_state.maps_Tf_inputvs.clear();
        for (int i = 0; i < n_meshes; i++)
            _map_state.maps_Tf_inputvs.push_back(ExternalProperty<FH, std::vector<SVH>> (_map_state.mesh_T));
    }

    // Compute vertex areas for input meshes
    if (_map_state.vertex_areas_input.size() == 0) // Only compute if not already done previously
    {
        for (int i = 0; i < n_meshes; i++)
            _map_state.vertex_areas_input.push_back(rel_vertex_areas(_map_state.meshes_input[i]));

        // If active vertices specified, scale active vertex areas to still reach total
        if ((_map_state.active_for_approx.size() == _map_state.meshes_input.size()))
        {
            for (int i = 0; i < n_meshes; i++)
            {
                double total_area = 0.0;
                double active_area = 0.0;
                for (auto vh : _map_state.meshes_input[i].vertices())
                {
                    total_area += _map_state.vertex_areas_input[0][vh];
                    if (_map_state.active_for_approx[i][vh])
                        active_area += _map_state.vertex_areas_input[0][vh];
                }
                for (auto vh : _map_state.meshes_input[i].vertices())
                {
                    if (_map_state.active_for_approx[i][vh])
                        _map_state.vertex_areas_input[0][vh] *= total_area / active_area;
                }
            }
        }
    }

    _callback_for_optim("init");

    // Alternate optimize_map() and remesh_T()
    double time_remeshing = 0.0;
    double time_optimize = 0.0;
    OptimizeMapStatus continuous_opt_status = normal;
    for (int iter = 0; iter < _settings.max_iterations; ++iter)
    {
        ISM_INFO("Starting iteration: " << iter);

        // Remeshing step
        bool remesh_changed_map_state = false;
        TinyAD::Timer timer_remeshing("remeshing");
        remesh_changed_map_state = remesh_T(_map_state, _settings, _csv_path, _callback_for_remesh);
        timer_remeshing.stop();
        time_remeshing += timer_remeshing.seconds();

        // Sanity checks
        ISM_ASSERT(_settings.w_approx <= 0.0 || sanity_check_all_vertices_assigned(_map_state));
        ISM_ASSERT(_settings.w_approx <= 0.0 || sanity_check_assignments(_map_state));

//        // Print stats
//        {
//            ISM_DEBUG_OUT("Vertices after remeshing: " << _map_state.mesh_T.n_vertices());
//            ISM_DEBUG_OUT("Faces after remeshing: " << _map_state.mesh_T.n_faces());
//            double min_el = INF_DOUBLE;
//            double max_el = 0.0;
//            double avg_el = 0.0;
//            std::vector<TriMesh> lifted_Ts = lifted_meshes_from_mapstate(_map_state);
//            for (int i = 0; i < n_meshes; i++)
//            {
//                for (auto e : _map_state.mesh_T.edges())
//                {
//                    const double el = lifted_Ts[i].calc_edge_length(e);
//                    min_el = std::min(min_el, el);
//                    max_el = std::max(min_el, el);
//                    avg_el += el;
//                }
//            }
//            avg_el /= n_meshes * _map_state.mesh_T.n_edges();
//            ISM_DEBUG_VAR(min_el);
//            ISM_DEBUG_VAR(max_el);
//            ISM_DEBUG_VAR(avg_el);
//        }

        // Check convergence after remeshing!
        if (!remesh_changed_map_state && continuous_opt_status == under_threshold)
        {
            ISM_HIGHLIGHT("Continous optimization and remeshing converged in iteration " << iter) ;
            break;
        }

        // No more progress?
        if (!remesh_changed_map_state && continuous_opt_status == line_search_converge)
        {
            ISM_HIGHLIGHT("Line search couldn't find improvement and remeshing converged in iteration " << iter);
            break;
        }

        // Landmarks satisfied?
        const double landmark_err = landmark_error(_map_state);
        ISM_DEBUG_VAR(landmark_err);
        if (!remesh_changed_map_state && landmark_error(_map_state) < _settings.landmark_threshold)
        {
            ISM_HIGHLIGHT("Landmark error under " << _settings.landmark_threshold << " in iteration " << iter);
            break;
        }

        // Continuous optimization
        TinyAD::Timer timer_optimize("optimize");
        continuous_opt_status = optimize_map(_map_state, _settings.optimize_per_remeshing_iters, _callback_for_optim, _csv_path, _settings);
        timer_optimize.stop();
        time_optimize += timer_optimize.seconds();

        // Sanity checks
        ISM_ASSERT(_settings.w_approx <= 0.0 || sanity_check_all_vertices_assigned(_map_state));
        ISM_ASSERT(_settings.w_approx <= 0.0 || sanity_check_assignments(_map_state));
    }
}

}

void optimize_with_remeshing(
        MapState& _map_state,
        const std::vector<ExternalProperty<VH, double>>& _mesh_tels,
        const AdaptiveTriangulationsSettings& _settings,
        const fs::path& _csv_path,
        std::function<void (const std::string&)> _callback_for_optim,
        std::function<void (const std::string&)> _callback_for_remesh)
{
    _map_state.tels_input = _mesh_tels;
    optimize_with_remeshing_impl(_map_state, _settings, _csv_path, _callback_for_optim, _callback_for_remesh);
}

void optimize_with_remeshing(
        MapState& _map_state,
        const AdaptiveTriangulationsSettings& _settings,
        const fs::path& _csv_path,
        std::function<void (const std::string&)> _callback_for_optim,
        std::function<void (const std::string&)> _callback_for_remesh)
{
    init_target_edge_lengths(_map_state, _settings);
    optimize_with_remeshing_impl(_map_state, _settings, _csv_path, _callback_for_optim, _callback_for_remesh);
}

}
