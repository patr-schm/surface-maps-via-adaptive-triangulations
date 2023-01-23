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
#include <SurfaceMaps/Utils/Genus.hh>

namespace SurfaceMaps
{

enum Geometry { Planar, Spherical, Hyperbolic };

Geometry which_geometry(const TriMesh& _mesh);

}
