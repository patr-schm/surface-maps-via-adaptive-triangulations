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
#include <SurfaceMaps/Utils/Geometries.hh>

namespace SurfaceMaps
{

bool flipped_or_degenerate(
        const Vec3d& _a,
        const Vec3d& _b,
        const Vec3d& _c,
        const Geometry& _geometry);

bool flipped_or_degenerate(
        const FH _fh,
        const TriMesh& _mesh,
        const ExternalProperty<VH, Vec3d>& _embedding,
        const Geometry& _geometry);

int count_flipped(
        const TriMesh& _mesh,
        const ExternalProperty<VH, Vec3d>& _embedding,
        const Geometry& _geometry);

}
