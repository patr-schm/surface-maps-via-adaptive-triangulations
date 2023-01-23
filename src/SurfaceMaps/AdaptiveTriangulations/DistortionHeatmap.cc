/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Dörte Pieper, Patrick Schmidt
 */

#include <SurfaceMaps/AdaptiveTriangulations/DistortionHeatmap.hh>

#include <SurfaceMaps/Utils/Helpers.hh>
#include <SurfaceMaps/Viewer/Colors.hh>
#include <SurfaceMaps/Viewer/HeatmapColors.hh>
#include <SurfaceMaps/AdaptiveTriangulations/Helpers.hh>

namespace SurfaceMaps
{

ExternalProperty<FH, Color> compute_distortion_heatmap(
        const TriMesh& _mesh_1,
        const TriMesh& _mesh_2,
        const float& _range_max)
{
    ISM_ASSERT_EQ(_mesh_1.n_faces(), _mesh_2.n_faces());

    ExternalProperty<FH, double> distortion(_mesh_1);
    float max_distortion = 0.0;

    for (auto fh : _mesh_1.faces())
    {
        Vec3d a_A, b_A, c_A, a_B, b_B, c_B;
        points(_mesh_1, fh, a_A, b_A, c_A);
        points(_mesh_2, fh, a_B, b_B, c_B);
        //ISM_DEBUG_OUT("fh " << fh);
        //local coordinate system per face
        Vec2d a_local_A, b_local_A, c_local_A, a_local_B, b_local_B, c_local_B;
        to_local_coordinates(a_A, b_A, c_A, a_local_A, b_local_A, c_local_A);
        to_local_coordinates(a_B, b_B, c_B, a_local_B, b_local_B, c_local_B);


        Eigen::Matrix2<double> M_A;
        M_A << b_local_A - a_local_A, c_local_A - a_local_A;
        Eigen::Matrix2<double> M_B;
        M_B << b_local_B - a_local_B, c_local_B - a_local_B;
        Eigen::Matrix2<double> J = M_B * M_A.inverse();
        Eigen::Matrix2<double> J_inv = M_A * M_B.inverse();

        distortion[fh] = J.squaredNorm() + J_inv.squaredNorm();
        if (distortion[fh] > max_distortion)
            max_distortion = distortion[fh];
    }

    // Convert (log) energy to colors
    const float range_min = 4.0f;
    const float range_max = _range_max;
    ISM_DEBUG_OUT("max distortion " << max_distortion);
    ExternalProperty<FH, Color> colors(_mesh_1);
    for (const FH fh_overlay :_mesh_1.faces())
        colors[fh_overlay] = log_color(distortion[fh_overlay], range_min, range_max, WHITE, MAGENTA);

    return colors;
}

}
