/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt, Dörte Pieper
 */
#include <SurfaceMaps/Init.hh>
#include <SurfaceMaps/Viewer/MeshView.hh>
#include <SurfaceMaps/AdaptiveTriangulations/InitSphereEmbeddings.hh>
#include <SurfaceMaps/AdaptiveTriangulations/OptimizeCoarseToFine.hh>
#include <SurfaceMaps/AdaptiveTriangulations/OptimizeWithRemeshing.hh>

const bool open_viewer = false; // or write screenshots

namespace SurfaceMaps
{

void run(
        const double w_approx)
{
    // Prepare output dir
    fs::path output_dir = OUTPUT_PATH / "approximation";
    fs::path screenshot_dir = output_dir / "screenshots";
    fs::create_directories(screenshot_dir);
    const std::string name = "w_approx_" + std::to_string(w_approx);

    const fs::path mesh_path_A = DATA_PATH / "meshes/hand/hand_plain.obj";
    const fs::path mesh_path_B = DATA_PATH / "meshes/hand/hand_vulcan.obj";
    const fs::path landmarks_path_A = DATA_PATH / "meshes/hand/hand_plain.pinned";
    const fs::path landmarks_path_B = DATA_PATH / "meshes/hand/hand_vulcan_0.pinned";
    const fs::path embedding_path_A = output_dir / "embedding_A.obj";
    const fs::path embedding_path_B = output_dir / "embedding_B.obj";

    // Init map
    MapState map_state;
    init_map(map_state, { mesh_path_A, mesh_path_B }, { landmarks_path_A, landmarks_path_B }, { embedding_path_A, embedding_path_B }, false);

    // Optimize map
    landmark_phase(map_state);
    release_landmarks(map_state);
    coarse_phase(map_state);

    AdaptiveTriangulationsSettings fine_settings = fine_phase_settings();
    fine_settings.w_approx = w_approx;
    optimize_with_remeshing(map_state, fine_settings, "");

    // Visualization
    std::vector<TriMesh> lifted_Ts = lifted_meshes_from_mapstate(map_state);

    // Screenshot T_A
    {
        auto cam_pos = glow::viewer::camera_transform(tg::pos3(1.732203f, 0.442215f, 3.076327f), tg::pos3(0.045052f, -0.113123f, -0.014236f));
        auto cam_config = gv::config(cam_pos);
        auto s = screenshot_config(open_viewer, screenshot_dir / (name + ".png"), tg::ivec2(1920, 1080), true);

        auto v = gv::view();
        auto style = gv::config(gv::sun_scale_factor(1.5));
        view_mesh(map_state.meshes_input[0], Color(1.0, 1.0, 1.0, 0.5));
        view_mesh(lifted_Ts[0], Color(1.0, 1.0, 1.0, 0.8));
        view_wireframe(lifted_Ts[0], MAGENTA, WidthScreen(0.7));
    }
}

}

int main()
{
    glow::glfw::GlfwContext ctx;
    using namespace SurfaceMaps;
    init_lib_surface_maps();

    for (double w_approx : { 0.0, 1.0 })
        run(w_approx);

    return 0;
}
