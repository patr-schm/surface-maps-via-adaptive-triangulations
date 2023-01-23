/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Dörte Pieper, Patrick Schmidt
 */

#include <SurfaceMaps/AdaptiveTriangulations/AdaptiveTargetEdgeLength.hh>

#include <SurfaceMaps/AdaptiveTriangulations/LiftToSurface.hh>
#include <SurfaceMaps/AdaptiveTriangulations/Helpers.hh>
#include <SurfaceMaps/Utils/PrincipalCurvature.hh>
#include <SurfaceMaps/Utils/IO.hh>

namespace SurfaceMaps
{

ExternalProperty<VH, double> target_edge_lengths(
        const TriMesh& _mesh,
        const double _approx_error,
        const double _min_edge_length,
        const double _max_edge_length)
{
    ISM_ASSERT_G(_approx_error, 0.0);
    ISM_ASSERT_L(_min_edge_length, _max_edge_length);

    const auto max_abs_curv = max_abs_curvatures(_mesh);

    ExternalProperty<VH, double> tels(_mesh);
    for (auto v : _mesh.vertices())
    {
        tels[v] = sqrt(6.0 * _approx_error / max_abs_curv[v] - 3.0 * sqr(_approx_error));
        tels[v] = std::clamp(tels[v], _min_edge_length, _max_edge_length);
        if(isnan(tels[v]))
            tels[v] = _min_edge_length;
    }

    return tels;
}

}
