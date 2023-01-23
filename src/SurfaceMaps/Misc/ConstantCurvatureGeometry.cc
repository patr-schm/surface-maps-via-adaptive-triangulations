/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt
 */

#include "ConstantCurvatureEmbedding.hh"

#include <ExactPredicates.h>

namespace SurfaceMaps
{

bool ccw_inclusive_3d(
        const Vec3d& a, const Vec3d& b, const Vec3d& c)
{
    Vec3d origin(0.0, 0.0, 0.0);
    return orient3d(a.data(), b.data(), c.data(), origin.data()) >= 0;
}

bool ccw_exclusive_3d(
        const Vec3d& a, const Vec3d& b, const Vec3d& c)
{
    Vec3d origin(0.0, 0.0, 0.0);
    return orient3d(a.data(), b.data(), c.data(), origin.data()) > 0;
}

bool in_triangle_inclusive_3d(
        const Vec3d& p,
        const Vec3d& a, const Vec3d& b, const Vec3d& c)
{
    return ccw_inclusive_3d(p, a, b) &&
           ccw_inclusive_3d(p, b, c) &&
           ccw_inclusive_3d(p, c, a);
}

bool in_triangle_exclusive_3d(
        const Vec3d& p,
        const Vec3d& a, const Vec3d& b, const Vec3d& c)
{
    return ccw_exclusive_3d(p, a, b) &&
           ccw_exclusive_3d(p, b, c) &&
           ccw_exclusive_3d(p, c, a);
}

}
