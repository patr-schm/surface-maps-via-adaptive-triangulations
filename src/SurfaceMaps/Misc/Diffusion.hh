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

template <typename T>
void diffuse_pointwise_field(
        const TriMesh& _mesh,
        ExternalProperty<VH, T>& _field,
        const double _t);

void diffuse_pointwise_field(
        const TriMesh& _mesh,
        ExternalProperty<VH, double>& _field,
        const double _t);
}

