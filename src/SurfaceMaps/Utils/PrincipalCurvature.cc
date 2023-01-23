/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt
 */

#include <SurfaceMaps/Utils/PrincipalCurvature.hh>

#include <SurfaceMaps/Utils/IO.hh>
#include <igl/principal_curvature.h>

namespace SurfaceMaps
{

ExternalProperty<VH, double> max_abs_curvatures(
        const TriMesh& _mesh,
        const int _k_ring)
{
    MatXd V;
    MatXi F;
    mesh_to_matrix(_mesh, V, F);

    MatXd max_dirs;
    MatXd min_dirs;
    VecXd max_curvs;
    VecXd min_curvs;
    igl::principal_curvature(V, F, max_dirs, min_dirs, max_curvs, min_curvs, _k_ring);

    ExternalProperty<VH, double> max_abs_curv(_mesh);
    max_abs_curv.as_eigen() = max_curvs.array().abs().max(min_curvs.array().abs());

    return max_abs_curv;
}

}
