/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt, Dörte Pieper
 */
#include <SurfaceMaps/Init.hh>
#include <SurfaceMaps/Utils/IO.hh>
#include <SurfaceMaps/Viewer/MeshView.hh>
#include <SurfaceMaps/Viewer/HeatmapColors.hh>
#include <SurfaceMaps/Utils/MeshNormalization.hh>

#include <SurfaceMaps/AdaptiveTriangulations/Helpers.hh>
#include <SurfaceMaps/AdaptiveTriangulations/Visualization.hh>
#include <SurfaceMaps/AdaptiveTriangulations/DistortionHeatmap.hh>
#include <SurfaceMaps/AdaptiveTriangulations/EvaluationMetrics.hh>
#include <SurfaceMaps/AdaptiveTriangulations/InitSphereEmbeddings.hh>
#include <SurfaceMaps/AdaptiveTriangulations/OptimizeCoarseToFine.hh>
#include <SurfaceMaps/AdaptiveTriangulations/OptimizeWithRemeshing.hh>

#include <TinyAD/Utils/Timer.hh>

const bool open_viewer = false; // or write screenshots

namespace SurfaceMaps
{

void run_pair(
        const double _approx_target,
        const fs::path _mesh_path_A,
        const fs::path _mesh_path_B,
        const fs::path _landmarks_path_A,
        const fs::path _landmarks_path_B,
        const std::string _name)
{
    ISM_DEBUG_OUT("Pair: " << _name);

    // Create output directory
    fs::path output_dir = OUTPUT_PATH / "resolution_T" / (_name);
    fs::path screenshot_dir = output_dir / ("screenshots_" + std::to_string(_approx_target));
    fs::create_directories(screenshot_dir);

    const fs::path embedding_path_A = output_dir / "embedding_A.obj";
    const fs::path embedding_path_B = output_dir / "embedding_B.obj";
    glow::SharedTexture2D texture = read_texture(DATA_PATH / "textures/checkerboard.png");

    // Init map
    TinyAD::Timer timer_init("Init (maybe cached)");
    MapState map_state;
    if (!init_map(map_state, { _mesh_path_A, _mesh_path_B }, { _landmarks_path_A, _landmarks_path_B }, { embedding_path_A, embedding_path_B }, false))
        return;
    timer_init.stop();

    rotate(map_state.meshes_input[0], -M_PI, Vec3d(0.0, 0.0, 1.0));
    rotate(map_state.meshes_input[1], -M_PI, Vec3d(0.0, 0.0, 1.0));

    // Visualization
    auto cam_pos = glow::viewer::camera_transform(tg::pos3(1.729787f, 0.130825f, 2.935483f), tg::pos3(1.241686f, 0.085696f, 1.963126f));
    double heatmap_max = 100.0;
    auto screenshot = [&] (const std::string& prefix)
    {
        TriMesh mesh_A_on_B = map_vertices_to_target(map_state, 0, 1);
        ExternalProperty<FH, Color> colors_AonB = compute_distortion_heatmap(map_state.meshes_input[0], mesh_A_on_B, heatmap_max);
        {
            // Mesh A on B
            auto cam_config = gv::config(cam_pos);
            auto s = screenshot_config(open_viewer, screenshot_dir / (prefix + "_AonB.png"), tg::ivec2(1920, 1080), true);
            auto v = gv::view();
            view_mesh(mesh_A_on_B, Color(1.0, 1.0, 1.0, 1.0));
            view_wireframe(mesh_A_on_B, MAGENTA, WidthScreen(0.7));
        }
        {
            // Heatmap
            auto cam_config = gv::config(cam_pos);
            auto s = screenshot_config(open_viewer, screenshot_dir / (prefix + "_AonB_heatmap.png"), tg::ivec2(1920, 1080), true);
            auto v = gv::view();
            auto style = default_style();
            gv::view(make_renderable(mesh_A_on_B, colors_AonB));
        }

        std::vector<TriMesh> lifted_Ts = lifted_meshes_from_mapstate(map_state);
        ExternalProperty<FH, Color> colors = compute_distortion_heatmap(lifted_Ts[0], lifted_Ts[1], heatmap_max);
        for (int i = 0; i < (int)map_state.meshes_input.size(); ++i)
        {
            auto sun = gv::config(gv::sun_scale_factor(1.5));
            {
                // Mesh T
                auto cam_config = gv::config(cam_pos);
                auto s = screenshot_config(open_viewer, screenshot_dir / (prefix + "_" + pad_integer(i, 2) + ".png"), tg::ivec2(1920, 1080), true);
                auto v = gv::view();
                view_mesh(lifted_Ts[i], Color(1.0, 1.0, 1.0, 1.0));
                view_wireframe(lifted_Ts[i], MAGENTA, WidthScreen(0.7));
                view_landmarks(lifted_Ts[i], map_state.landmarks_T, WidthScreen(10.0));
            }

            {
                // Texture
                auto cam_config = gv::config(cam_pos);
                auto s = screenshot_config(open_viewer, screenshot_dir / (prefix + "_" + pad_integer(i, 2) + "_textured.png"), tg::ivec2(1920, 1080), true);
                auto v = gv::view();
                view_texture_frontal_projection_input(map_state, 0, i, 2, 1.5, texture);
                view_landmarks(lifted_Ts[i], map_state.landmarks_T, WidthScreen(10.0));
            }
            {
                // Heatmap
                auto cam_config = gv::config(cam_pos);
                auto s = screenshot_config(open_viewer, screenshot_dir / (prefix + "_" + pad_integer(i, 2) + "_heatmap.png"), tg::ivec2(1920, 1080), true);
                auto v = gv::view();
                auto style = default_style();
                gv::view(make_renderable(lifted_Ts[i], colors));
            }
        }
    };

    // Algorithm schedule
    screenshot("00_init");
    TinyAD::Timer timer_landmark("Landmark Phase");
    landmark_phase(map_state);
    timer_landmark.stop();
    screenshot("01_after_landmark");

    release_landmarks(map_state);

    double diag_A = bounding_box_diagonal(map_state.meshes_input[0]);
    double diag_B = bounding_box_diagonal(map_state.meshes_input[1]);
    double final_approx_error = _approx_target * fmin(diag_A, diag_B);
    ISM_DEBUG_OUT("final approx error: " << final_approx_error);
    ISM_DEBUG_OUT("coarse phase error: " << fmax(final_approx_error, 0.01))

    TinyAD::Timer timer_coarse("Coarse Phase");
    coarse_phase(map_state, fmax(final_approx_error, 0.01));
    screenshot("02_after_coarse");
    timer_coarse.stop();

    TinyAD::Timer timer_fine("Fine Phase");

    if (final_approx_error <= 0.01)
    {
        fine_phase(map_state, final_approx_error);
        screenshot("03_after_fine");
    }
    timer_fine.stop();

    const double runtime_seconds = timer_init.seconds() + timer_landmark.seconds() + timer_coarse.seconds() + timer_fine.seconds();
    ISM_INFO("Total run time: " << runtime_seconds << " seconds");

    // Save lifted T meshes
    {
        const std::vector<TriMesh> liftedTs = lifted_meshes_from_mapstate(map_state);
        write_mesh(liftedTs[0], output_dir / ("T_A_" + std::to_string(_approx_target) + ".obj"));
        write_mesh(liftedTs[1], output_dir / ("T_B_" + std::to_string(_approx_target) + ".obj"));
    }

    // Evaluate ours
    {
        std::vector<TriMesh> lifted_Ts = lifted_meshes_from_mapstate(map_state);
        evaluate_meshes(map_state.meshes_input[0], lifted_Ts[0], map_state.meshes_input[1], lifted_Ts[1], runtime_seconds, output_dir / ("eval_" + std::to_string(_approx_target) + ".txt"));

        // Write distortion to csv
        write_mapping_distortion_normalized(lifted_Ts[0], lifted_Ts[1], output_dir / ("distortion_" + std::to_string(_approx_target) + ".csv"));
    }

    // Distortion between input A and AonB
    {
        TriMesh mesh_A_on_B = map_vertices_to_target(map_state, 0, 1);
        // Write distortion to csv
        write_mapping_distortion_normalized(map_state.meshes_input[0], mesh_A_on_B, output_dir / ("distortion_AonB_" + std::to_string(_approx_target) + ".csv"));
    }
}

}

int main()
{
    glow::glfw::GlfwContext ctx;
    using namespace SurfaceMaps;
    init_lib_surface_maps();

    const fs::path mesh_path_A = DATA_PATH / "meshes/yang2020/cat_horse/input0.obj";
    const fs::path mesh_path_B = DATA_PATH / "meshes/yang2020/cat_horse/input1.obj";
    const fs::path landmarks_path_A = DATA_PATH / "meshes/yang2020/cat_horse/landmarks0.txt";
    const fs::path landmarks_path_B = DATA_PATH / "meshes/yang2020/cat_horse/landmarks1.txt";

    run_pair(0.04, mesh_path_A, mesh_path_B, landmarks_path_A, landmarks_path_B, "cat_horse");
    run_pair(0.024, mesh_path_A, mesh_path_B, landmarks_path_A, landmarks_path_B, "cat_horse");
    run_pair(0.016, mesh_path_A, mesh_path_B, landmarks_path_A, landmarks_path_B, "cat_horse");
    run_pair(0.008, mesh_path_A, mesh_path_B, landmarks_path_A, landmarks_path_B, "cat_horse");
    run_pair(0.004, mesh_path_A, mesh_path_B, landmarks_path_A, landmarks_path_B, "cat_horse");
    run_pair(0.002, mesh_path_A, mesh_path_B, landmarks_path_A, landmarks_path_B, "cat_horse");
    run_pair(0.001, mesh_path_A, mesh_path_B, landmarks_path_A, landmarks_path_B, "cat_horse");

    return 0;
}
