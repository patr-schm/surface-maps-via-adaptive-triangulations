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
#include <SurfaceMaps/AdaptiveTriangulations/OptimizeWithRemeshing.hh>
#include <SurfaceMaps/Utils/MeshNormalization.hh>
#include <SurfaceMaps/Misc/Diffusion.hh>

#include <TinyAD/Utils/Timer.hh>

const bool open_viewer = false; // or write screenshots

namespace SurfaceMaps
{

void run_pair(
        const int _idx_A,
        const int _idx_B,
        const std::string _scalarfield_name_A,
        const std::string _scalarfield_name_B)
{
    ISM_DEBUG_OUT("Pair: A " << _idx_A << " B " << _idx_B);
    std::string name = std::to_string(_idx_A) + "_" + std::to_string(_idx_B);

    // Create output directory
    fs::path output_dir = OUTPUT_PATH / "user_sizing" / name;
    fs::path screenshot_dir = output_dir / "screenshots";
    fs::create_directories(screenshot_dir);

    const fs::path mesh_path_A = DATA_PATH / "meshes/shrec" / (std::to_string(_idx_A) + ".off");
    const fs::path mesh_path_B = DATA_PATH / "meshes/shrec" / (std::to_string(_idx_B) + ".off");
    const fs::path landmarks_path_A = DATA_PATH / "landmarks/shrec_user_sizing" / (std::to_string(_idx_A) + ".vts");
    const fs::path landmarks_path_B = DATA_PATH / "landmarks/shrec_user_sizing" / (std::to_string(_idx_B) + ".vts");
    const fs::path embedding_path_A = output_dir / "embedding_A.obj";
    const fs::path embedding_path_B = output_dir / "embedding_B.obj";

    MapState map_state;
    if (!init_map(map_state, { mesh_path_A, mesh_path_B }, { landmarks_path_A, landmarks_path_B }, { embedding_path_A, embedding_path_B }, true))
        return;

    rotate(map_state.meshes_input[0], -M_PI/18.0, Vec3d(1.0, 0.0, 0.0));
    rotate(map_state.meshes_input[1], -M_PI/8.0, Vec3d(1.0, 0.0, 0.0));

    // Visualization
    int screenshot_idx = 0;
    auto screenshot = [&] (const std::string& comment)
    {
        const fs::path screenshot_path = screenshot_dir / (pad_integer(screenshot_idx) + "_" + comment + ".png");

        auto s = screenshot_config(open_viewer, screenshot_path, tg::ivec2(1920, 1080));
        view_lifted_with_landmarks_in_meshes(map_state);
        ++screenshot_idx;
    };

    screenshot("init");

    // Optimize (landmark and coarse phase)
    landmark_phase(map_state);
    screenshot("after_landmarks");
    coarse_phase(map_state, 0.002);
    screenshot("after_coarse");

    // Load scalar field
    map_state.tels_input[0] = read_property<VH>(DATA_PATH / "fields" / name / _scalarfield_name_A, map_state.meshes_input[0], true);
    map_state.tels_input[1] = read_property<VH>(DATA_PATH / "fields" / name / _scalarfield_name_B, map_state.meshes_input[1], true);
    diffuse_pointwise_field(map_state.meshes_input[0], map_state.tels_input[0], 1e-4);

    // Optimize (fine phase)
    AdaptiveTriangulationsSettings fine_settings = fine_phase_settings();
    fine_settings.w_barrier = 5e-6;
    optimize_with_remeshing(map_state, map_state.tels_input, fine_settings, "");
    screenshot("after_fine");

    // Screenshots
    std::vector<TriMesh> lifted_Ts = lifted_meshes_from_mapstate(map_state);
    auto cam_pos = glow::viewer::camera_transform(tg::pos3(-0.835719f, 0.556923f, 1.950492f), tg::pos3(-0.000482f, 0.036468f, -0.014201f));
    {
        auto s = screenshot_config(open_viewer, screenshot_dir/ "scalarfield.png", cam_pos,  tg::ivec2(1920, 1080), true);
        {
            view_scalar_field(map_state.meshes_input[0], map_state.tels_input[0], RED, BLUE, true);
        }
    }
    {
        {
            auto v = gv::view();
            auto s = screenshot_config(open_viewer, screenshot_dir/ "wireframeMeshesA.png", cam_pos,  tg::ivec2(1920, 1080), true);
            view_mesh(lifted_Ts[0]);
            view_wireframe(lifted_Ts[0], MAGENTA);
        }
        {
            auto cam_pos2 = glow::viewer::camera_transform(tg::pos3(0.635719f, 0.556923f, 1.950492f), tg::pos3(-0.000482f, 0.036468f, -0.014201f));
            auto v = gv::view();
            auto s = screenshot_config(open_viewer, screenshot_dir/ "wireframeMeshesB.png", cam_pos2,  tg::ivec2(1920, 1080), true);
            view_mesh(lifted_Ts[1]);
            view_wireframe(lifted_Ts[1], MAGENTA);
        }
    }

    // Write result meshes
    for (int i = 0; i < (int)map_state.meshes_input.size(); ++i)
        write_mesh(lifted_Ts[i], output_dir / ("T_" + pad_integer(i, 2) + ".obj"));
}

}

int main()
{
    glow::glfw::GlfwContext ctx;
    using namespace SurfaceMaps;
    init_lib_surface_maps();

    run_pair(302, 314, "field_" + std::to_string(302) + "_A", "field_" + std::to_string(314) + "_B");

    return 0;
}
