/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt, Janis Born
 */

#include "LandmarkEditor.hh"

#include <imgui/imgui.h>
#include <SurfaceMaps/Utils/IO.hh>
#include <SurfaceMaps/Utils/Helpers.hh>
#include <SurfaceMaps/Utils/Dijkstra.hh>
#include <SurfaceMaps/Utils/MeshNormalization.hh>
#include <SurfaceMaps/Viewer/MeshView.hh>
#include <SurfaceMaps/Viewer/Picking.hh>
#include <SurfaceMaps/Viewer/ColorGenerator.hh>

#include <glow-extras/viewer/canvas.hh>

namespace SurfaceMaps
{

namespace
{

void update_points_and_lines(
        LandmarkEditor& _editor)
{
    _editor.draw.clear();

    {
        ColorGenerator cg;
        for (const auto& vh : _editor.landmarks)
            _editor.draw.point(_editor.mesh.point(vh), cg.generate_next_color(), WidthScreen(_editor.point_size_px, true));
    }

    {
        ColorGenerator cg;
        for (const auto& segment : _editor.paths)
        {
            const auto c = cg.generate_next_color();
            for (auto& e : segment)
                _editor.draw.line(e, _editor.mesh, c);
        }
    }
}

bool vertex_picked(
        const LandmarkEditor& _editor,
        const VH& _vh)
{
    for (auto v : _editor.landmarks)
    {
        if (v == _vh)
            return true;
    }

    return false;
}

int vertex_idx(
        const LandmarkEditor& _editor,
        const VH& _vh)
{
    for (int idx = 0; idx < (int)_editor.landmarks.size(); ++idx)
    {
        if (_editor.landmarks[idx] == _vh)
            return idx;
    }

    return -1;
}

bool edge_picked(
        const LandmarkEditor& _editor,
        const EH& _eh)
{
    for (auto& segment : _editor.paths)
    {
        for (auto e : segment)
        {
            if (e == _eh)
                return true;
        }
    }

    return false;
}

void add_shortest_path(
        LandmarkEditor& _editor,
        const VH _vh_from,
        const VH _vh_to)
{
    ExternalProperty<EH, bool> blocked_edges(_editor.mesh, false);
    for (auto& segment : _editor.paths)
    {
        for (auto e : segment)
            blocked_edges[e] = true;
    }

    auto path = find_primal_path(_vh_from, _vh_to, _editor.mesh, blocked_edges, ExternalProperty<VH, bool>(_editor.mesh, false));

    _editor.paths.push_back({});
    for (auto h : path.hehs)
        _editor.paths.back().push_back(_editor.mesh.edge_handle(h));
}

void render_gui(
        LandmarkEditor& _editor)
{
    // Create the editor GUI
    ImGui::Begin("Landmark and Path Editor");

    // Display section
    ImGui::Separator();
    ImGui::Text("Display");
    ImGui::Separator();
    if (ImGui::SliderFloat("Point size", &_editor.point_size_px, 1.0, 16.0, "%.1f"))
        update_points_and_lines(_editor);
    ImGui::Checkbox("Show vertex IDs", &_editor.show_landmark_ids);

    // Picking section
    ImGui::Separator();
    ImGui::Text("Picking");
    ImGui::Separator();

    // Toggle pick modes
    if (ImGui::Checkbox("Pick vertices", &_editor.pick_vertices) && _editor.pick_vertices)
        _editor.pick_edges = false;
    if (ImGui::Checkbox("Pick edges", &_editor.pick_edges) && _editor.pick_edges)
        _editor.pick_vertices = false;

    ImGui::Checkbox("Add shortest paths", &_editor.add_shortest_paths);

    if (ImGui::Button("Close path loop"))
    {
        if (_editor.landmarks.size() >= 2)
            add_shortest_path(_editor, _editor.landmarks.back(), _editor.landmarks.front());

        update_points_and_lines(_editor);
    }

    // Pick vertices (middle-mouse shortcut)
    if (ImGui::IsMouseClicked(ImGuiMouseButton_Middle))
    {
        const VH vh = pick_vertex(_editor.mesh);
        if (vh.is_valid())
        {
            // Add path from previous to this landmark
            if (_editor.add_shortest_paths && !_editor.landmarks.empty())
                add_shortest_path(_editor, _editor.landmarks.back(), vh);

            // Add new landmark
            _editor.landmarks.push_back(vh);

            update_points_and_lines(_editor);
        }
    }

    // Delete vertices (right-mouse shortcut)
    if (ImGui::IsMouseClicked(ImGuiMouseButton_Right))
    {
        if (_editor.landmarks.size() > 0)
        {
            _editor.landmarks.pop_back();
            update_points_and_lines(_editor);
        }
    }

    // Pick vertices or edges
    if (_editor.pick_vertices)
    {
        if (ImGui::IsMouseClicked(ImGuiMouseButton_Left))
        {
            const VH vh = pick_vertex(_editor.mesh);
            if (vh.is_valid())
            {
                if (vertex_picked(_editor, vh))
                {
                    // Start dragging existing landmark
                    _editor.dragging_vertex_idx = vertex_idx(_editor, vh);
                }
                else
                {
                    // Add new landmark and start dragging
                    _editor.landmarks.push_back(vh);
                    _editor.dragging_vertex_idx = _editor.landmarks.size() - 1;
                }

                update_points_and_lines(_editor);
            }
        }

        if (ImGui::IsMouseDragging(ImGuiMouseButton_Left))
        {
            if (_editor.dragging_vertex_idx >= 0)
            {
                const VH vh = pick_vertex(_editor.mesh);
                if (vh.is_valid())
                {
                    // Relocate dragged vertex
                    _editor.landmarks[_editor.dragging_vertex_idx] = vh;
                }
                else
                {
                    // Stop dragging
                    _editor.dragging_vertex_idx = -1;
                }

                update_points_and_lines(_editor);
            }
        }

        if (ImGui::IsMouseReleased(ImGuiMouseButton_Left))
        {
            if (_editor.dragging_vertex_idx >= 0)
            {
                // Stop dragging
                _editor.dragging_vertex_idx = -1;

                update_points_and_lines(_editor);
            }
        }
    }

    if (_editor.pick_edges)
    {
        if (ImGui::IsMouseClicked(ImGuiMouseButton_Left))
        {
            auto eh = pick_edge(_editor.mesh);
            if (eh.is_valid() && ! edge_picked(_editor, eh))
            {
                _editor.paths.push_back({ eh });
                update_points_and_lines(_editor);
            }
        }
    }

    // Clearing section
    ImGui::Separator();
    ImGui::Text("Clear");
    ImGui::Separator();

    if (ImGui::Button("Remove last vertex"))
    {
        if (_editor.landmarks.size() > 0)
        {
            _editor.landmarks.pop_back();
            update_points_and_lines(_editor);
        }
    }
    if (_editor.landmarks.size() > 0)
    {
        ImGui::SameLine();
        ImGui::Text("(ID %d)", _editor.landmarks.back().idx());
    }

    if (ImGui::Button("Remove last edge segment"))
    {
        if (_editor.paths.size() > 0)
        {
            _editor.paths.pop_back();
            update_points_and_lines(_editor);
        }
    }
    if (_editor.paths.size() > 0 && _editor.paths.back().size() > 0)
    {
        ImGui::SameLine();
        ImGui::Text("(ID %d)", _editor.paths.back().back().idx());
    }

    if (ImGui::Button("Clear vertices"))
    {
        _editor.landmarks.clear();
        update_points_and_lines(_editor);
    }
    if (ImGui::Button("Clear edges"))
    {
        _editor.paths.clear();
        update_points_and_lines(_editor);
    }

    // Load/Save section
    ImGui::Separator();
    ImGui::Text("Save");
    ImGui::Separator();

    ImGui::InputInt("Extra init landmarks start at", &_editor.init_landmarks_begin_idx);

    ImGui::Checkbox("Overwrite files", &_editor.overwrite);

    if (ImGui::Button("Save landmarks"))
    {
        std::vector<LandmarkType> landmark_types;
        for (int i = 0; i < (int)_editor.landmarks.size(); ++i)
        {
            if (i < _editor.init_landmarks_begin_idx)
                landmark_types.push_back(LandmarkType::Keep);
            else
                landmark_types.push_back(LandmarkType::Init);
        }

        write_landmarks(_editor.landmarks, _editor.path_landmarks, landmark_types, _editor.overwrite);
    }

    if (ImGui::Button("Save paths"))
    {
        std::vector<EH> ehs;
        for (auto& segment : _editor.paths)
        {
            for (auto e : segment)
                ehs.push_back(e);
        }
        write_paths(ehs, _editor.path_paths, _editor.overwrite);
    }

    ImGui::End();
}

}

LandmarkEditor::LandmarkEditor(
        const fs::path& _path_mesh,
        const fs::path& _path_landmarks,
        const fs::path& _path_paths,
        const fs::path& _path_texture)
    : path_mesh(_path_mesh),
      path_landmarks(_path_landmarks),
      path_paths(_path_paths)
{
    // Load mesh
    mesh = read_mesh(_path_mesh);
    normalize_surface_area(mesh);

    // Load landmarks
    if (path_landmarks.empty())
        path_landmarks = path_mesh.parent_path() / (path_mesh.stem().string() + ".pinned");
    landmarks = read_landmarks(path_landmarks, LandmarkType::Init, false);

    // Load paths
    if (path_paths.empty())
        path_paths = path_landmarks.parent_path() / (path_landmarks.stem().string() + ".paths");
    paths.push_back(read_paths(path_paths, false));

    // Load texture
    if (_path_texture.empty())
        r_mesh = make_renderable(mesh);
    else
    {
        glow::SharedTexture2D texture = glow::Texture2D::createFromFile(_path_texture, glow::ColorSpace::sRGB);
        r_mesh = make_renderable(mesh, texture);
    }

    update_points_and_lines(*this);
}

void LandmarkEditor::view()
{
    // Disable mouse camera control when picking.
    // This has to go before the first gv::view() call.
    GV_SCOPED_CONFIG(gv::no_left_mouse_control(pick_vertices || pick_edges));

    auto v = view_mesh(r_mesh);
    render_gui(*this);
    draw.view();

    if (show_landmark_ids)
    {
        auto c = gv::canvas();
        for (size_t i = 0; i < landmarks.size(); ++i)
        {
            const VH lm = landmarks[i];
            const Vec3d pos = mesh.point(lm);
            const tg::pos3 pos_tg = tg::pos3(pos[0], pos[1], pos[2]);
            c.add_label(pos_tg, std::to_string(i));
        }
    }
}

}
