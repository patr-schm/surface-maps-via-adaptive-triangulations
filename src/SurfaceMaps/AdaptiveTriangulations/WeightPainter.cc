/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Dörte Pieper, Patrick Schmidt
 */

#include <SurfaceMaps/AdaptiveTriangulations/WeightPainter.hh>
#include <SurfaceMaps/Utils/IO.hh>
#include <SurfaceMaps/Viewer/Picking.hh>
#include <SurfaceMaps/Viewer/ColorGenerator.hh>
#include <imgui/imgui.h>
#include <glow-extras/viewer/canvas.hh>
#include <glow-extras/viewer/experimental.hh>
#include <SurfaceMaps/Utils/IO.hh>

#include <SurfaceMaps/Viewer/MeshView.hh>

#include <SurfaceMaps/AdaptiveTriangulations/Visualization.hh>
#include <SurfaceMaps/AdaptiveTriangulations/Helpers.hh>

namespace SurfaceMaps
{
FH pick_face(const TriMesh& _mesh)
{
    using namespace gv::experimental;

    // Don't pick if UI captures the mouse click
    if (ImGui::GetIO().WantCaptureMouse)
        return FH(-1);

    // Get mouse position
    auto p_world = interactive_get_position(interactive_get_mouse_position());

    if (!p_world.has_value() || _mesh.n_faces() == 0)
        return FH(-1);

    // Return closest edge midpoint
    const Vec3d p(p_world.value()[0], p_world.value()[1], p_world.value()[2]);
    return _mesh.faces().argmin([&] (FH fh) { return (_mesh.calc_face_centroid(fh) - p).squaredNorm(); });

}


std::vector<VH> pick_vertices_in_radius(
        const TriMesh& _mesh,
        const float& radius)
{
    std::vector<VH> picked_vertices;


    // Don't pick if UI captures the mouse click
    if (ImGui::GetIO().WantCaptureMouse)
        return picked_vertices;

    // Get mouse position
    auto p_world = gv::experimental::interactive_get_position(gv::experimental::interactive_get_mouse_position());
    if (!p_world.has_value() || _mesh.n_vertices() == 0)
        return picked_vertices;

    const auto p = to_vec3<double>(p_world.value());

    FH picked_face = pick_face(_mesh);
    Vec3d face_normal = _mesh.calc_face_normal(picked_face);

    for (auto vh : _mesh.vertices())
    {
        // Collect all in radius
        if ((_mesh.point(vh) - p).norm() < radius)
        {
            Vec3d normal_vh;
            _mesh.calc_vertex_normal_correct(vh, normal_vh);
            // Only paint on front
            if(normal_vh.dot(face_normal) > 0)
                picked_vertices.push_back(vh);
        }
    }
    return picked_vertices;

}

void view_weight_painter(
        const TriMesh& _mesh,
        ExternalProperty<VH, double>& _tels,
        const fs::path& _target_dir_scalar)
{
    if (_tels.empty() )
        _tels = ExternalProperty<VH, double> (_mesh, 0.0);

    // gui state
    double target_scalar_value = 0.0;
    float radius = 0.01;
    int softness = 0;
    std::vector<tg::segment3> circle;

    int max_softness = 100;
    // View
    std::vector<VH> change_vertices;

    gv::interactive([&] (auto)
    {
        ImGui::Begin("Options");
        ImGui::InputDouble("Value", &target_scalar_value);
        ImGui::SliderFloat("Radius", &radius, 0.01, 1.0);
        ImGui::SliderInt("Softness", &softness, 0, max_softness);
        if(!_target_dir_scalar.empty())
        {
            if (ImGui::Button("Save Sizing Field"))
            {
                write_property(_tels, _target_dir_scalar, true);
            }
        }

        ImGui::End();

        auto const mouse_pos = gv::experimental::interactive_get_mouse_position();
        auto const pick = gv::experimental::interactive_get_position(mouse_pos);
        if (pick.has_value() && !ImGui::GetIO().WantCaptureMouse)
        {
            auto const pick_pos = tg::pos3(pick.value());

            // Pick face and calculate normal
            Vec3d face_normal, pa, pb, pc;
            FH picked_face = pick_face(_mesh);
            face_normal = _mesh.calc_face_normal(picked_face);
            points(_mesh, picked_face, pa, pb, pc);

            // Tangent for draw circle
            Vec3d face_tangent1 = (pa-pb).normalized();
            Vec3d face_tangent2 = face_tangent1.cross(face_normal).normalized();

            auto const n_segs = 32;
            circle.clear();
            for (auto i = 0; i < n_segs; ++i)
            {
                auto const a0 = float(i) * 360_deg / float(n_segs);
                auto const a1 = float(i + 1) * 360_deg / float(n_segs);
                auto const s0 = tg::sin(a0);
                auto const c0 = tg::cos(a0);
                auto const s1 = tg::sin(a1);
                auto const c1 = tg::cos(a1);
                auto const p0 = tg::pos3(pick_pos) + tg::vec3(face_tangent1) * c0 * radius + tg::vec3(face_tangent2) * s0 * radius;
                auto const p1 = tg::pos3(pick_pos) + tg::vec3(face_tangent1) * c1 * radius + tg::vec3(face_tangent2) * s1 * radius;
                circle.push_back({p0, p1});
            }
        }

        // Pick on the mesh
        if (ImGui::IsMouseDown(ImGuiMouseButton_Middle))
        {
            change_vertices = pick_vertices_in_radius(_mesh, radius);
        }

        // View
        auto g = gv::grid();

        {
            auto sub_view = gv::view();
            // Draw annotations
            ColorGenerator colors(2);
            auto c = gv::canvas();

            for (auto vh : change_vertices)
            {
                _tels[vh] = (target_scalar_value * (max_softness - softness) + _tels[vh] * softness)/double(max_softness);
                change_vertices.erase(change_vertices.begin());
            }
            gv::view(gv::lines(circle).camera_facing(), tg::color3::black, gv::maybe_empty);

            view_caption("A");
            view_scalar_field(_mesh, _tels, RED, BLUE);

            sub_view.configure(gv::no_shadow);
        }
    });

}

}
