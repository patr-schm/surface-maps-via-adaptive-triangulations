/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Joe Jakobi, Patrick Schmidt
 */
#pragma once

#include <SurfaceMaps/Types.hh>
#include <SurfaceMaps/MultiRes/MultiResSphereEmbeddingSettings.hh>

namespace SurfaceMaps
{

struct MultiResSphereEmbedding
{
    MultiResSphereEmbedding(
            TriMesh& _mesh,
            int _n_vertices_total,
            OpenMesh::VPropHandleT<Vec3d>& _embedding,
            ProgressiveMesh& _pm,
            const std::vector<uint>& _independent_set_sizes)
        : mesh(_mesh),
          n_vertices_total(_n_vertices_total),
          ph_embedding(_embedding),
          pm(_pm),
          independent_set_sizes(_independent_set_sizes)
    {
        rec_it = _pm.log.crbegin();
        current_independent_set = _independent_set_sizes.size() - 1;
    };

    TriMesh& mesh;
    int n_vertices_total;
    const OpenMesh::VPropHandleT<Vec3d> ph_embedding;
    ProgressiveMesh pm;
    const std::vector<uint> independent_set_sizes;
    uint current_independent_set;
    std::vector<ProgressiveMesh::Record>::const_reverse_iterator rec_it;
    double mesh_area = 0.0;
};

ExternalProperty<VH, Vec3d> multi_res_sphere_embedding(
        const TriMesh& _mesh,
        const MultiResSphereEmbeddingSettings& _settings = MultiResSphereEmbeddingSettings(),
        std::function<void(const TriMesh&, const ExternalProperty<VH, Vec3d>&)> _callback = [] (const TriMesh&, const ExternalProperty<VH, Vec3d>&) {});
}
