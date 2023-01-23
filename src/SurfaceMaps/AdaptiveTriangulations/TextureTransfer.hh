/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Janis Born, Patrick Schmidt
 */
#pragma once

#include <SurfaceMaps/AdaptiveTriangulations/MapState.hh>
#include <glow-extras/viewer/view.hh>

namespace SurfaceMaps
{

gv::SharedRenderable make_renderable_smooth(
        const TriMesh& _mesh,
        const glow::SharedTexture2D& _texture);

TriMesh transfer_texture_and_subdiv(
        const MapState& _map_state,
        const int _idx_from,
        const int _idx_to,
        const int _n_subdiv = 2);

}
