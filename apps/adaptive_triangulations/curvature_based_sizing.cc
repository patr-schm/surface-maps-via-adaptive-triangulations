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
#include <SurfaceMaps/AdaptiveTriangulations/Visualization.hh>
#include <SurfaceMaps/AdaptiveTriangulations/InitSphereEmbeddings.hh>
#include <SurfaceMaps/AdaptiveTriangulations/OptimizeCoarseToFine.hh>

const bool open_viewer = false; // or write screenshots

int main()
{
    glow::glfw::GlfwContext ctx;
    using namespace SurfaceMaps;
    init_lib_surface_maps();

    // Prepare output dir
    fs::path output_dir = OUTPUT_PATH / "curvature_based_sizing";
    fs::path screenshot_dir = output_dir / "screenshots";
    fs::create_directories(screenshot_dir);

    const fs::path mesh_path_A = DATA_PATH / "meshes/yang2020/rabbit12-15/input0.obj";
    const fs::path mesh_path_B = DATA_PATH / "meshes/yang2020/rabbit12-15/input1.obj";
    const fs::path landmarks_path_A = DATA_PATH / "meshes/yang2020/rabbit12-15/landmarks0.txt";
    const fs::path landmarks_path_B = DATA_PATH / "meshes/yang2020/rabbit12-15/landmarks1.txt";
    const fs::path embedding_path_A = output_dir / "embedding_A.obj";
    const fs::path embedding_path_B = output_dir / "embedding_B.obj";

    // Init map
    MapState map_state;
    init_map(map_state, { mesh_path_A, mesh_path_B }, { landmarks_path_A, landmarks_path_B }, { embedding_path_A, embedding_path_B }, false);

    rotate(map_state.meshes_input[0], -M_PI / 1.0, Vec3d(1.0, 0.0, 0.0));
    rotate(map_state.meshes_input[1], -M_PI / 1.0, Vec3d(1.0, 0.0, 0.0));

    auto cam_pos = glow::viewer::camera_transform(tg::pos3(1.458124f, 0.324316f, -1.990489f), tg::pos3(-0.015362f, -0.025248f, 0.001074f));
    auto screenshot = [&] (const std::string& prefix)
    {
        std::vector<TriMesh> lifted_Ts = lifted_meshes_from_mapstate(map_state);
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
                // Scalar Field target edge length
                auto cam_config = gv::config(cam_pos);
                auto s = screenshot_config(open_viewer, screenshot_dir / (prefix + "_" + pad_integer(i, 2) + "_scalarfield.png"), tg::ivec2(1920, 1080), true);
                {
                    TriMesh _mesh = map_state.meshes_input[i];
                    auto _field = map_state.tels_input[i];
                    auto _color_from = RED;
                    auto _color_to = BLUE;

                    auto style = default_style();
                    pm::Mesh m;
                    auto pos = to_polymesh(_mesh, m);

                    auto v_colors = m.vertices().make_attribute<tg::color3>();
                    const double min = 0.007;
                    const double max = 0.05;

                    const double range = max - min;

                    tg::color3 tg_color_from(_color_from[0], _color_from[1], _color_from[2]);
                    tg::color3 tg_color_to(_color_to[0], _color_to[1], _color_to[2]);

                    for (auto v : m.vertices())
                    {
                        if (range <= 0.0)
                        {
                            v_colors[v] = tg_color_from;
                            continue;
                        }

                        const double val = _field[VH(v.idx.value)];
                        v_colors[v] = tg::color3(log_color(val, min, max, _color_from, _color_to));
                    }

                    gv::view(pos, v_colors);
                }
            }
            {
                // Texture
                glow::SharedTexture2D texture = read_texture(DATA_PATH / "textures/checkerboard.png");
                auto cam_config = gv::config(cam_pos);
                auto s = screenshot_config(open_viewer, screenshot_dir / (prefix + "_" + pad_integer(i, 2) + "_textured.png"), tg::ivec2(1920, 1080), true);
                auto v = gv::view();
                view_texture_frontal_projection_input(map_state, 0, i, 2, 1.5, texture);
                view_landmarks(lifted_Ts[i], map_state.landmarks_T, WidthScreen(10.0));
            }
        }
    };

    // Optimize map
    landmark_phase(map_state);
    release_landmarks(map_state);
    coarse_phase(map_state);
    screenshot("01_after_coarse");
    fine_phase(map_state);
    screenshot("02_after_fine");

    // Write result meshes
    std::vector<TriMesh> lifted_Ts = lifted_meshes_from_mapstate(map_state);
    for (int i = 0; i < (int)map_state.meshes_input.size(); ++i)
        write_mesh(lifted_Ts[i], output_dir / ("T_" + pad_integer(i, 2) + ".obj"));

    return 0;
}
