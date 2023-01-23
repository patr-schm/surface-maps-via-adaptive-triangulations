/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt
 */
#pragma once

#include <SurfaceMaps/Types.hh>

namespace SurfaceMaps
{

ExternalProperty<VH, double> max_abs_curvatures(
        const TriMesh& _mesh,
        const int _k_ring = 5);

}
