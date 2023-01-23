/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt, Janis Born
 */

# include "Picking.hh"

#include <imgui/imgui.h>
#include <SurfaceMaps/Utils/Out.hh>
#include <SurfaceMaps/Utils/Helpers.hh>
#include <glow-extras/viewer/experimental.hh>

namespace SurfaceMaps
{

Vec3d to_vec3(const tg::pos3& _p)
{
    return Vec3d(_p.x, _p.y, _p.z);
}

VH pick_vertex(
        const TriMesh& _mesh)
{
    using namespace gv::experimental;

    // Don't pick if UI captures the mouse click
    if (ImGui::GetIO().WantCaptureMouse)
        return VH(-1);

    // Get mouse position
    auto p_world = interactive_get_position(interactive_get_mouse_position());

    if (!p_world.has_value() || _mesh.n_vertices() == 0)
        return VH(-1);

    // Return closest vertex
    const auto p = to_vec3<double>(p_world.value());
    return _mesh.vertices().argmin([&] (VH vh) { return (_mesh.point(vh) - p).squaredNorm(); });
}

EH pick_edge(
        const TriMesh& _mesh)
{
    using namespace gv::experimental;

    // Don't pick if UI captures the mouse click
    if (ImGui::GetIO().WantCaptureMouse)
        return EH(-1);

     // Get mouse position
     auto p_world = interactive_get_position(interactive_get_mouse_position());

    if (!p_world.has_value() || _mesh.n_edges() == 0)
        return EH(-1);

    // Return closest edge midpoint
    const auto p = to_vec3<double>(p_world.value());
    return _mesh.edges().argmin([&] (EH eh) { return (_mesh.calc_edge_midpoint(eh) - p).squaredNorm(); });
}

}
