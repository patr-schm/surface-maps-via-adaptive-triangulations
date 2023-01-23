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
#include <SurfaceMaps/Utils/Filesystem.hh>
#include <SurfaceMaps/AdaptiveTriangulations/MapState.hh>
#include <SurfaceMaps/MultiRes/MultiResSphereEmbedding.hh>

namespace SurfaceMaps
{

bool init_map(
        MapState& _map_state,
        const std::vector<fs::path>& _mesh_paths,
        const std::vector<fs::path>& _landmark_paths,
        const std::vector<fs::path>& _embedding_paths,
        const bool _align_meshes = true);

void compute_input_sphere_embeddings(
        const TriMesh& _mesh_A,
        const TriMesh& _mesh_B,
        const std::vector<VH>& _landmarks_A,
        const std::vector<VH>& _landmarks_B,
        ExternalProperty<VH, Vec3d>& _embedding_A,
        ExternalProperty<VH, Vec3d>& _embedding_B,
        const MultiResSphereEmbeddingSettings& _settings = MultiResSphereEmbeddingSettings());


}
