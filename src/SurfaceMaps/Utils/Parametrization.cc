/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt
 */

#include "Parametrization.hh"

#include <SurfaceMaps/Utils/Helpers.hh>
#include <ExactPredicates.h>

namespace SurfaceMaps
{

bool ccw_exact_exclusive_2d(
        const Vec2d& a, const Vec2d& b, const Vec2d& c)
{
    return orient2d(a.data(), b.data(), c.data()) > 0.0;
}

bool flipped_or_degenerate(
        const Vec2d& _a,
        const Vec2d& _b,
        const Vec2d& _c)
{
    return !ccw_exact_exclusive_2d(_a, _b, _c);
}

bool flipped_or_degenerate(
        const FH _fh,
        const TriMesh& _mesh,
        const Parametrization& _param)
{
    Vec2d a, b, c;
    points(_mesh, _param, _fh, a, b, c);

    return flipped_or_degenerate(a, b, c);
}

}
