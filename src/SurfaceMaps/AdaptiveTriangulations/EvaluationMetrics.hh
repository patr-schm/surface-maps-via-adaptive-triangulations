/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Dörte Pieper, Patrick Schmidt
 */
#pragma once

#include <SurfaceMaps/Types.hh>
#include <SurfaceMaps/Utils/Filesystem.hh>
#include <SurfaceMaps/AdaptiveTriangulations/MapState.hh>

namespace SurfaceMaps
{

TriMesh map_vertices_to_target(
        const MapState& _map_state,
        const int _source_mesh_idx,
        const int _target_mesh_idx);

///one sided hausdorff distance points mesh_1 to surface mesh_2, samples on vertices, edges and faces
double compute_hausdorff_distance(
        const TriMesh& _mesh_1,
        const TriMesh& _mesh_2,
        const int& _samples_per_face = 1);

/// two sided hausdorff distances between mesh_orig to lifted mesh_approx
double hausdorff_distance_for_trimeshes(
        const TriMesh& _mesh_orig,
        const TriMesh& _mesh_approx);

void mapping_distortion_normalized(
        const TriMesh& _mesh_TA,
        const TriMesh& _mesh_TB,
        ExternalProperty<FH, double>& _areas_TA,
        ExternalProperty<FH, double>& _areas_TB,
        ExternalProperty<FH, double>& _sing_vals_min,
        ExternalProperty<FH, double>& _sing_vals_max);

void write_mapping_distortion_normalized(
        const TriMesh& _mesh_TA,
        const TriMesh& _mesh_TB,
        const fs::path& _csv_path);

/// isometric distortion via singular values of jacobian
void metric_mapping_distortion(
        const TriMesh &_mesh_1,
        const TriMesh &_mesh_2,
        double & _distortion);

/// evaluate meshes with different metrics
void evaluate_meshes(
        const TriMesh& _mesh_A,
        const TriMesh& _mesh_A_approx,
        const TriMesh& _mesh_B,
        const TriMesh& _mesh_B_approx,
        const double _runtime_seconds,
        const fs::path& _out_path="",
        const bool& _silent = false);

}
