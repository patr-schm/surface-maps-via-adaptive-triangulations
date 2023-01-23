/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Dörte Pieper, Patrick Schmidt
 */
#pragma once

#include <SurfaceMaps/Types.hh>
#include <SurfaceMaps/AdaptiveTriangulations/MapState.hh>

namespace SurfaceMaps
{

double max_dist_verts_T(
        const MapState& _map_state);

double n_verts_over_bound(
        const MapState& _map_state,
        const double& _bound);

void satisfy_approx_bound_face_splits(
        MapState& _map_state,
        const double &_bound,
        const int& _max_iters = 20);

void local_refine_tel(
        MapState& _map_state,
        const double& _bound,
        const double& _factor,
        const bool _split);

}
