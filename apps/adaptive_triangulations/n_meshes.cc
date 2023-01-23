/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt, Dörte Pieper
 */
#include <SurfaceMaps/Init.hh>
#include <SurfaceMaps/Utils/Genus.hh>
#include <SurfaceMaps/Utils/IO.hh>
#include <SurfaceMaps/Viewer/MeshView.hh>
#include <SurfaceMaps/Utils/MeshNormalization.hh>
#include <SurfaceMaps/AdaptiveTriangulations/Visualization.hh>
#include <SurfaceMaps/AdaptiveTriangulations/InitSphereEmbeddings.hh>
#include <SurfaceMaps/AdaptiveTriangulations/OptimizeCoarseToFine.hh>
#include <SurfaceMaps/AdaptiveTriangulations/OptimizeWithRemeshing.hh>
#include <TinyAD/Utils/Timer.hh>

const bool open_viewer = false; // or write screenshots

namespace SurfaceMaps
{

void run_n_meshes(
        const std::vector<int> _indices,
        const std::string _subdir_name = "test")
{
    std::string dir_name = "n_meshes/" + _subdir_name;
    std::string file_name = _subdir_name;

    // Create output directory
    fs::path output_dir = OUTPUT_PATH / dir_name;
    fs::path screenshot_dir = output_dir / "screenshots";
    fs::create_directories(screenshot_dir);

    // Collect paths to meshes, landmarks, and sphere embeddings
    std::vector<fs::path> mesh_paths;
    std::vector<fs::path> landmark_paths;
    std::vector<fs::path> embedding_paths;
    for (int i = 0; i < (int)_indices.size(); i++)
    {
        mesh_paths.push_back(DATA_PATH / "meshes/shrec" / (std::to_string(_indices[i]) + ".off"));
        landmark_paths.push_back(DATA_PATH / "meshes/shrec" / (std::to_string(_indices[i]) + ".vts"));
        embedding_paths.push_back(output_dir / ("embedding_" + std::to_string(_indices[i]) +".obj"));
    }
    glow::SharedTexture2D texture = read_texture(DATA_PATH / "textures/checkerboard.png");

    const bool clear_cached = false;
    if (clear_cached)
    {
        for (auto& path : embedding_paths)
            fs::remove_all(path);
    }

    // Load meshes and init map
    MapState map_state;
    TinyAD::Timer timer_init("Init (maybe cached)");
    init_map(map_state, mesh_paths, landmark_paths, embedding_paths, false);
    timer_init.stop();

    map_state.set_distortion_pairs(DistortionPairs::Star);

    // Rotate meshes
    const int ref_align = 0;
    for (int i = 0; i < (int)_indices.size(); i++)
    {
        std::vector<VH> feet_landmarks_ref = {map_state.landmarks_input[ref_align][5], map_state.landmarks_input[ref_align][7], map_state.landmarks_input[ref_align][9], map_state.landmarks_input[ref_align][11]};
        std::vector<VH> feet_landmarks_i = {map_state.landmarks_input[i][5], map_state.landmarks_input[i][7], map_state.landmarks_input[i][9], map_state.landmarks_input[i][11]};
        align_rigid(feet_landmarks_ref, feet_landmarks_i, map_state.meshes_input[ref_align], map_state.meshes_input[i]);
    }

    // Visualization
    auto cam_pos = glow::viewer::camera_transform(tg::pos3(0.091298f, 0.816960f, 4.112880f), tg::pos3(0.111494f, 0.671299f, 3.258824f));
    auto screenshot = [&] (const std::string& prefix)
    {
        std::vector<TriMesh> lifted_Ts = lifted_meshes_from_mapstate(map_state);
        for (int i = 0; i < (int)map_state.meshes_input.size(); ++i)
        {
            // View T mesh
            {
                auto cam_config = gv::config(cam_pos);
                auto s = screenshot_config(open_viewer, screenshot_dir / (prefix + "_" + pad_integer(i, 2) + ".png"), tg::ivec2(1920, 1080), true);

                auto v = gv::view();
                auto style = gv::config(gv::sun_scale_factor(2));
                view_mesh(map_state.meshes_input[i], Color(1.0, 1.0, 1.0, 0.2));
                view_mesh(lifted_Ts[i], Color(1.0, 1.0, 1.0, 0.9));
                view_wireframe(lifted_Ts[i], MAGENTA, WidthScreen(0.6));
                view_landmarks(lifted_Ts[i], map_state.landmarks_T, WidthScreen(10.0));
            }

            // View texture
            {
                auto cam_config = gv::config(cam_pos);
                auto s = screenshot_config(open_viewer, screenshot_dir / (prefix + "_" + pad_integer(i, 2) + "_textured.png"), tg::ivec2(1920, 1080), true);

                auto v = gv::view();
                auto style = gv::config(gv::sun_scale_factor(2));
                view_texture_frontal_projection_input(map_state, 0, i, 2, 0.6, texture);
                view_landmarks(lifted_Ts[i], map_state.landmarks_T, WidthScreen(10.0));
            }
        }
    };

    // Optimize
    screenshot("00_init");

    TinyAD::Timer timer_landmark("Landmark Phase");
    landmark_phase(map_state);
    timer_landmark.stop();

    screenshot("01_after_landmark");

    release_landmarks(map_state);

    TinyAD::Timer timer_coarse("Coarse Phase");
    AdaptiveTriangulationsSettings coarse_settings = coarse_phase_settings();
    coarse_settings.w_approx = 4.0;
    optimize_with_remeshing(map_state, coarse_settings);
    timer_coarse.stop();

    screenshot("02_after_coarse");

    TinyAD::Timer timer_fine("Fine Phase");
    AdaptiveTriangulationsSettings fine_settings = fine_phase_settings(0.01);
    fine_settings.max_iterations = 30;
    fine_settings.w_approx = 4.0;
    optimize_with_remeshing(map_state, fine_settings);
    timer_fine.stop();

    screenshot("03_after_fine");

    // Write result meshes
    std::vector<TriMesh> lifted_Ts = lifted_meshes_from_mapstate(map_state);
    for (int i = 0; i < (int)map_state.meshes_input.size(); ++i)
        write_mesh(lifted_Ts[i], output_dir / ("T_" + pad_integer(i, 2) + ".obj"));

    ISM_INFO("Total run time: " << timer_init.seconds() + timer_landmark.seconds() + timer_coarse.seconds() + timer_fine.seconds() << " seconds");
}

}

int main()
{
    glow::glfw::GlfwContext ctx;
    using namespace SurfaceMaps;
    init_lib_surface_maps();

    run_n_meshes({ 381, 382, 383, 384, 385, 387, 388, 389, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400 }, "fourleg");

    return 0;
}
