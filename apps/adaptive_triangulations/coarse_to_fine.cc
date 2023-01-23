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
#include <SurfaceMaps/AdaptiveTriangulations/Visualization.hh>
#include <SurfaceMaps/AdaptiveTriangulations/InitSphereEmbeddings.hh>
#include <SurfaceMaps/AdaptiveTriangulations/OptimizeCoarseToFine.hh>
#include <TinyAD/Utils/Timer.hh>

const bool open_viewer = false; // or write screenshots

int main()
{
    glow::glfw::GlfwContext ctx;
    using namespace SurfaceMaps;
    init_lib_surface_maps();

    // Prepare output dir
    fs::path output_dir = OUTPUT_PATH / "coarse_to_fine";
    fs::path screenshot_dir = output_dir / "screenshots";
    fs::create_directories(screenshot_dir);

    const fs::path mesh_path_A = DATA_PATH / "meshes/shrec/61.off";
    const fs::path mesh_path_B = DATA_PATH / "meshes/shrec/248.off";
    const fs::path landmarks_path_A = DATA_PATH / "landmarks/shrec_plane_bird/61.pinned";
    const fs::path landmarks_path_B = DATA_PATH / "landmarks/shrec_plane_bird/248.pinned";
    const fs::path embedding_path_A = output_dir / "embedding_A.obj";
    const fs::path embedding_path_B = output_dir / "embedding_B.obj";
    glow::SharedTexture2D texture = read_texture(DATA_PATH / "textures/checkerboard.png");

    const bool clear_cached = false;
    if (clear_cached)
    {
        fs::remove_all(embedding_path_A);
        fs::remove_all(embedding_path_B);
    }

    // Load meshes and init map
    TinyAD::Timer timer_init("Init (maybe cached)");
    MapState map_state;
    init_map(map_state, { mesh_path_A, mesh_path_B }, { landmarks_path_A, landmarks_path_B }, { embedding_path_A, embedding_path_B }, false);
    timer_init.stop();

    // Rotate meshes
    rotate(map_state.meshes_input[0], -M_PI / 2.0, Vec3d(1.0, 0.0, 0.0));
    rotate(map_state.meshes_input[1], -M_PI / 2.0, Vec3d(0.0, 1.0, 0.0));

    // Visualization
    auto cam_pos = glow::viewer::camera_transform(tg::pos3(-1.400362f, 2.122391f, 0.559440f), tg::pos3(-0.868983f, 1.343911f, 0.346269f));
    auto screenshots = [&] (const std::string& prefix)
    {
        // View T meshes
        std::vector<TriMesh> lifted_Ts = lifted_meshes_from_mapstate(map_state);
        {
            auto cam_config = gv::config(cam_pos);
            auto s = screenshot_config(open_viewer, screenshot_dir / (prefix + "_A.png"), tg::ivec2(1920, 1080), true);

            auto v = gv::view();
            auto style = gv::config(gv::sun_scale_factor(1.5));
            view_mesh(map_state.meshes_input[0], Color(0.85, 0.85, 0.85, 0.5));
            view_mesh(lifted_Ts[0], Color(1.0, 1.0, 1.0, 0.8));
            view_wireframe(lifted_Ts[0], MAGENTA, WidthScreen(0.7));
            view_landmarks(lifted_Ts[0], map_state.landmarks_T, WidthScreen(25.0));
        }
        {
            auto cam_config = gv::config(cam_pos);
            auto s = screenshot_config(open_viewer, screenshot_dir / (prefix + "_B.png"), tg::ivec2(1920, 1080), true);

            auto v = gv::view();
            auto style = gv::config(gv::sun_scale_factor(1.5));
            view_mesh(map_state.meshes_input[1], Color(0.85, 0.85, 0.85, 0.5));
            view_mesh(lifted_Ts[1], Color(1.0, 1.0, 1.0, 0.8));
            view_wireframe(lifted_Ts[1], MAGENTA, WidthScreen(0.7));
            view_landmarks(lifted_Ts[1], map_state.landmarks_T, WidthScreen(20.0));
        }

        // View T wireframes without shadow (top layer)
        {
            auto cam_config = gv::config(cam_pos);
            auto s = screenshot_config(open_viewer, screenshot_dir / (prefix + "_A_wire.png"), tg::ivec2(1920, 1080), true);

            auto v = gv::view();
            auto style = gv::config(gv::no_shadow);
            view_mesh(lifted_Ts[0], Color(1.0, 1.0, 1.0, 0.8));
            view_wireframe(lifted_Ts[0], MAGENTA, WidthScreen(0.7));
            view_landmarks(lifted_Ts[0], map_state.landmarks_T, WidthScreen(25.0));
        }
        {
            auto cam_config = gv::config(cam_pos);
            auto s = screenshot_config(open_viewer, screenshot_dir / (prefix + "_B_wire.png"), tg::ivec2(1920, 1080), true);

            auto v = gv::view();
            auto style = gv::config(gv::no_shadow);
            view_mesh(lifted_Ts[1], Color(1.0, 1.0, 1.0, 0.8));
            view_wireframe(lifted_Ts[1], MAGENTA, WidthScreen(0.7));
            view_landmarks(lifted_Ts[1], map_state.landmarks_T, WidthScreen(20.0));
        }

        // View textures
        {
            auto cam_config = gv::config(cam_pos);
            auto s = screenshot_config(open_viewer, screenshot_dir / (prefix + "_A_textured.png"), tg::ivec2(1920, 1080), true);

            auto v = gv::view();
            auto style = gv::config(gv::sun_scale_factor(1.5));
            view_texture_frontal_projection_input(map_state, 0, 0, 1, 1.0, texture);
            view_landmarks(lifted_Ts[0], map_state.landmarks_T, WidthScreen(25.0));
        }
        {
            auto cam_config = gv::config(cam_pos);
            auto s = screenshot_config(open_viewer, screenshot_dir / (prefix + "_B_textured.png"), tg::ivec2(1920, 1080), true);

            auto v = gv::view();
            auto style = gv::config(gv::sun_scale_factor(1.5));
            view_texture_frontal_projection_input(map_state, 0, 1, 1, 1.0, texture);
            view_landmarks(lifted_Ts[1], map_state.landmarks_T, WidthScreen(20.0));
        }
    };

    // Optimize map
    screenshots("00_init");

    TinyAD::Timer timer_landmark("Landmark Phase");
    landmark_phase(map_state);
    timer_landmark.stop();
    screenshots("01_after_landmark");

    TinyAD::Timer timer_coarse("Coarse Phase");
    coarse_phase(map_state);
    timer_coarse.stop();
    screenshots("02_after_coarse");

    TinyAD::Timer timer_fine("Fine Phase");
    fine_phase(map_state);
    timer_fine.stop();
    screenshots("03_after_fine");

    // Write result meshes
    std::vector<TriMesh> lifted_Ts = lifted_meshes_from_mapstate(map_state);
    for (int i = 0; i < (int)map_state.meshes_input.size(); ++i)
        write_mesh(lifted_Ts[i], output_dir / ("T_" + pad_integer(i, 2) + ".obj"));

    ISM_INFO("Total run time: " << timer_init.seconds() + timer_landmark.seconds() + timer_coarse.seconds() + timer_fine.seconds() << " seconds");

    return 0;
}
