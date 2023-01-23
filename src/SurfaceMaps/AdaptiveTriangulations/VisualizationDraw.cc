/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Dörte Pieper, Patrick Schmidt
 */

#include <SurfaceMaps/AdaptiveTriangulations/VisualizationDraw.hh>

#include <SurfaceMaps/Viewer/ColorGenerator.hh>
#include <SurfaceMaps/Viewer/Picking.hh>
#include <SurfaceMaps/Utils/IO.hh>
#include <glow-extras/viewer/canvas.hh>
#include <imgui/imgui.h>

#include <SurfaceMaps/Viewer/MeshView.hh>
#include <SurfaceMaps/AdaptiveTriangulations/Visualization.hh>
#include <SurfaceMaps/AdaptiveTriangulations/Helpers.hh>

namespace SurfaceMaps
{

void draw_arrow(
        const TriMesh& _mesh,
        const VH _vh,
        const Color& _color,
        gv::canvas_t& c)
{
    const double arrow_length = 0.07;
    const double arrow_width = 0.01;
    auto to = _mesh.point(_vh);
    auto from = to + arrow_length * _mesh.calc_normal(_vh);
    c.add_arrow(tg::pos3(from), tg::pos3(to), arrow_width, tg::color3(_color));
}

void draw_line(
        const TriMesh& _mesh,
        const std::vector<VH>& _vhs,
        const Color& _color,
        gv::canvas_t& c)
{
    for (int i = 0; i < (int)_vhs.size() - 1; ++i)
    {
        auto from = _mesh.point(_vhs[i]);
        auto to = _mesh.point(_vhs[i+1]);
        c.set_line_width_world(0.01);
        c.add_line(tg::pos3(from), tg::pos3(to), tg::color3(_color));
    }
}

struct Annotation
{
    std::vector<VH> vhs; // on mesh_T
};

void view_lifted_draw(
        MapState& _map_state)
{
    const int n_meshes = _map_state.meshes_input.size();
    const std::vector<TriMesh> lifted_Ts = lifted_meshes_from_mapstate(_map_state);

    // View
    std::vector<Annotation> annotations;
    bool dragging = false;
    gv::interactive([&] (auto)
    {
        // Pick only in first view
        VH vh_pick = pick_vertex(lifted_Ts[0]);

        // Start arrow/line
        if (ImGui::IsMouseDown(ImGuiMouseButton_Middle))
        {
            if (vh_pick.is_valid())
            {
                if (!dragging)
                    annotations.push_back(Annotation {{ vh_pick }});
                else if (annotations.back().vhs.back() != vh_pick)
                    annotations.back().vhs.push_back(vh_pick);
            }

            dragging = true;
        }

        // Finish arrow/line
        if (ImGui::IsMouseReleased(ImGuiMouseButton_Middle))
            dragging = false;

        // Delete latest annotation
        if (ImGui::IsMouseClicked(ImGuiMouseButton_Right) && !annotations.empty())
            annotations.erase(annotations.end() - 1);

        // View
        auto g = gv::grid();

        for (int i = 0; i < n_meshes; ++i)
        {
            auto sub_view = gv::view();
            auto c = gv::canvas();

            // Draw annotations
            ColorGenerator colors(2);
            for (const auto& a : annotations)
            {
                if (a.vhs.size() == 1)
                {
                    draw_arrow(lifted_Ts[i], a.vhs.front(), colors.generate_next_color(), c);
                }
                else
                {
                    draw_line(lifted_Ts[i], a.vhs, colors.generate_next_color(), c);
                }
            }

            view_mesh(lifted_Ts[i]);
            view_wireframe(lifted_Ts[i], MAGENTA);
        }
    });
}

}
