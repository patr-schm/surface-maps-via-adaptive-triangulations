/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt
 */

#include "ConstantCurvatureEmbedding.hh"

#include <SurfaceMaps/Utils/Helpers.hh>
#include <SurfaceMaps/Misc/ConstantCurvatureGeometry.hh>
#include <ExactPredicates.h>
#include <queue>

namespace SurfaceMaps
{

bool flipped_or_degenerate(
        const Vec3d& _a,
        const Vec3d& _b,
        const Vec3d& _c,
        const Geometry& _geometry)
{
    ISM_ASSERT(_geometry == Spherical || _geometry == Planar || _geometry == Hyperbolic);

    // Perform both inexact and checks for consistency
    Mat3d M;
    M << _a, _b, _c;
    if (M.determinant() <= 0.0)
        return true;

    if (!ccw_exclusive_3d(_a, _b, _c))
        return true;

    return false;
}

bool flipped_or_degenerate(
        const FH _fh,
        const TriMesh& _mesh,
        const ExternalProperty<VH, Vec3d>& _embedding,
        const Geometry& _geometry)
{
    VH vh_a, vh_b, vh_c;
    handles(_mesh, _fh, vh_a, vh_b, vh_c);

    return flipped_or_degenerate(_embedding[vh_a], _embedding[vh_b], _embedding[vh_c], _geometry);
}

int count_flipped(
        const TriMesh& _mesh,
        const ExternalProperty<VH, Vec3d>& _embedding,
        const Geometry& _geometry)
{
    int n = 0.0;
    for (auto f : _mesh.faces())
    {
        if (flipped_or_degenerate(f, _mesh, _embedding, _geometry))
            ++n;
    }

    return n;
}

}
