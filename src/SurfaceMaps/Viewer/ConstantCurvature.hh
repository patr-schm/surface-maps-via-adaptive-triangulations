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
#include <SurfaceMaps/Utils/Geometries.hh>
#include <SurfaceMaps/Viewer/IDraw.hh>
#include <SurfaceMaps/Viewer/Colors.hh>

#include <glow-extras/glfw/GlfwContext.hh>
#include <glow-extras/viewer/view.hh>

namespace SurfaceMaps
{

gv::detail::raii_view_closer view_ccm_embedding(
        const TriMesh& _mesh,
        const ExternalProperty<VH, Vec3d>& _embedding,
        const Geometry& _geometry,
        const Color& _color = BLUE,
        const bool _view_manifold = true,
        const int _n_samples = 100,
        const DrawStyle& _style = default_line_style);

}
