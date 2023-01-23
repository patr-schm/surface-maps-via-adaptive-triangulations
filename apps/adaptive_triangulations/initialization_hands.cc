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

#include <SurfaceMaps/AdaptiveTriangulations/Helpers.hh>
#include <SurfaceMaps/AdaptiveTriangulations/Visualization.hh>
#include <SurfaceMaps/AdaptiveTriangulations/DistortionHeatmap.hh>
#include <SurfaceMaps/AdaptiveTriangulations/InitSphereEmbeddings.hh>
#include <SurfaceMaps/AdaptiveTriangulations/OptimizeCoarseToFine.hh>
#include <SurfaceMaps/AdaptiveTriangulations/OptimizeWithRemeshing.hh>
#include <SurfaceMaps/AdaptiveTriangulations/LiftToSurface.hh>

const bool open_viewer = false; // or write screenshots

namespace SurfaceMaps
{

void run_pair(
        const int _handB_num)
{
    // Create output directory
    fs::path output_dir = OUTPUT_PATH / "initialization_hands" / ("hand_" + std::to_string(_handB_num));
    fs::path screenshot_dir = output_dir / ("screenshots");
    fs::create_directories(screenshot_dir);

    const fs::path _mesh_path_A = DATA_PATH / "meshes/hand/hand_plain.obj";
    const fs::path _mesh_path_B = DATA_PATH / "meshes/hand/hand_vulcan.obj";
    const fs::path _landmarks_path_A = DATA_PATH / "meshes/hand/hand_plain.pinned";
    const fs::path _landmarks_path_B = DATA_PATH / ("meshes/hand/hand_vulcan_" + std:: to_string(_handB_num) + ".pinned");
    const fs::path embedding_path_A = output_dir / "embedding_A.obj";
    const fs::path embedding_path_B = output_dir / "embedding_B.obj";
    glow::SharedTexture2D texture = read_texture(DATA_PATH / "meshes/hand/hand_plain_texture.png");

    // Init map
    MapState map_state;
    if (!init_map(map_state, { _mesh_path_A, _mesh_path_B }, { _landmarks_path_A, _landmarks_path_B }, { embedding_path_A, embedding_path_B }, false))
        return;

    // Visualization
    auto cam_pos = glow::viewer::camera_transform(tg::pos3(-0.155488f, 0.903014f, 1.873404f), tg::pos3(-0.128843f, 0.788787f, 1.659081f));
    double heatmap_max = 100.0;

    auto screenshot = [&] (const std::string& prefix)
    {
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
                view_mesh(map_state.meshes_input[i], Color(1.0, 1.0, 1.0, 0.5));
                view_mesh(lifted_Ts[i], Color(1.0, 1.0, 1.0, 1.0));
                view_wireframe(lifted_Ts[i], MAGENTA, WidthScreen(0.7));
                view_landmarks(lifted_Ts[i], map_state.landmarks_T, WidthScreen(10.0));
            }

            {
                // Texture
                auto cam_config = gv::config(cam_pos);
                auto s = screenshot_config(open_viewer, screenshot_dir / (prefix + "_" + pad_integer(i, 2) + "_textured.png"), tg::ivec2(1920, 1080), true);
                auto v = gv::view();
                view_texture_halfedgetexcoords_input(map_state, 0, i, texture);
                view_landmarks(lifted_Ts[i], map_state.landmarks_T, WidthScreen(10.0));
            }
        }
    };

    ISM_ASSERT(map_state.meshes_input[0].has_halfedge_texcoords2D());

    // Algorithm schedule
    screenshot("00_init");
    landmark_phase(map_state);
    screenshot("01_after_landmark");

    release_landmarks(map_state);

    AdaptiveTriangulationsSettings coarse_settings = coarse_phase_settings();
    coarse_settings.w_map = 10.0;
    coarse_settings.set_adaptive_tel_params(0.015, 0.001, 100.0);
    coarse_settings.max_iterations = 400;
    optimize_with_remeshing(map_state, coarse_settings, "");
    screenshot("02_after_coarse");

    AdaptiveTriangulationsSettings coarse_settings2 = coarse_phase_settings();
    coarse_settings2.max_iterations = 200;
    optimize_with_remeshing(map_state, coarse_settings2, "");
    screenshot("03_after_coarse");

    AdaptiveTriangulationsSettings fine_settings = fine_phase_settings(0.0008);
    fine_settings.max_iterations = 200;
    optimize_with_remeshing(map_state, fine_settings, "");
    screenshot("04_after_fine");

    // Save lifted T meshes
    {
        const std::vector<TriMesh> liftedTs = lifted_meshes_from_mapstate(map_state);
        write_mesh(liftedTs[0], output_dir / ("T_A.obj"));
        write_mesh(liftedTs[1], output_dir / ("T_B.obj"));
    }
}

}

int main()
{
    glow::glfw::GlfwContext ctx;
    using namespace SurfaceMaps;
    init_lib_surface_maps();

    run_pair(1);
    run_pair(3);
    run_pair(5);
    run_pair(7);
    run_pair(8);

    return 0;
}
