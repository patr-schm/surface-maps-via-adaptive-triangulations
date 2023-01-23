/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt, Janis Born
 */

#pragma once

#include <polymesh/formats.hh>
#include <glow-extras/viewer/view.hh>
#include <SurfaceMaps/Types.hh>

namespace SurfaceMaps
{

/// Pick closest vertex in mesh
VH pick_vertex(
        const TriMesh& _mesh);

/// Pick closest edge midpoint in mesh
EH pick_edge(
        const TriMesh& _mesh);

}
