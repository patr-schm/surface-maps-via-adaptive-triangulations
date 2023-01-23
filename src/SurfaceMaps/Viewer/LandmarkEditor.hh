/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt, Janis Born
 */

#pragma once

#include <SurfaceMaps/Types.hh>
#include <SurfaceMaps/Viewer/MeshView.hh>
#include <SurfaceMaps/Viewer/GlowDraw.hh>

namespace SurfaceMaps
{

struct LandmarkEditor
{
    LandmarkEditor(
            const fs::path& _path_mesh,
            const fs::path& _path_landmarks = fs::path(),
            const fs::path& _path_paths = fs::path(),
            const fs::path& _path_texture = fs::path());

    void view();

    bool show_landmark_ids = false;
    float point_size_px = 8.0f;

    bool pick_vertices = false;
    bool pick_edges = false;
    bool add_shortest_paths = false;

    int dragging_vertex_idx = -1;

    fs::path path_mesh;
    fs::path path_landmarks;
    fs::path path_paths;

    TriMesh mesh;
    std::vector<VH> landmarks;
    std::vector<std::vector<EH>> paths;

    int init_landmarks_begin_idx = std::numeric_limits<int>::max();
    bool overwrite = false;

    gv::SharedRenderable r_mesh;
    GlowDraw draw;
};

}
