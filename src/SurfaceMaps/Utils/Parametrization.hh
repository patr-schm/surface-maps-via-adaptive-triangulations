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

/**
 * Performs exact check for single triangle.
 * True if non-positive area.
 */
bool flipped_or_degenerate(
        const Vec2d& _a,
        const Vec2d& _b,
        const Vec2d& _c);

/**
 * Performs exact check for single triangle.
 * True if non-positive area.
 */
bool flipped_or_degenerate(
        const FH _fh,
        const TriMesh& _mesh,
        const Parametrization& _param);

}
