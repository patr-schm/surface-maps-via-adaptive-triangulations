/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Dörte Pieper, Patrick Schmidt
 */

#include <SurfaceMaps/AdaptiveTriangulations/MapState.hh>

#include <SurfaceMaps/AdaptiveTriangulations/Helpers.hh>
#include <SurfaceMaps/AdaptiveTriangulations/LiftToSurface.hh>

namespace SurfaceMaps
{

void MapState::compute_bsp_trees()
{
    ISM_ASSERT_EQ(embeddings_input.size(), meshes_input.size());

    meshes_embeddings_input.resize(meshes_input.size());
    bsp_embeddings_input.resize(meshes_input.size());
    for (int i = 0; i < (int)meshes_input.size(); i++)
    {
        meshes_embeddings_input[i] = embedding_to_mesh(meshes_input[i], embeddings_input[i]);
        bsp_embeddings_input[i] = (BSPTree(meshes_embeddings_input[i]));
    }
}

void MapState::set_sphere_embeddings(
        const std::vector<ExternalProperty<VH, Vec3d>>& _embeddings)
{
    ISM_ASSERT_EQ(_embeddings.size(), meshes_input.size());

    embeddings_input.resize(meshes_input.size());
    for (int i = 0; i < (int)meshes_input.size(); i++)
        embeddings_input[i] = _embeddings[i];

    compute_bsp_trees();
}

void MapState::set_distortion_pairs(
        const DistortionPairs& _pairs_mode)
{
    pairs_map_distortion.clear();

    const int n_meshes = meshes_input.size();
    if (_pairs_mode == DistortionPairs::None)
    {

    }
    else if (_pairs_mode == DistortionPairs::All)
    {
        for (int i = 0; i < n_meshes; ++i)
        {
            for (int j = i + 1; j < n_meshes; ++j)
                pairs_map_distortion.push_back({ i, j });
        }
    }
    else if (_pairs_mode == DistortionPairs::Sequence)
    {
        for (int i = 0; i < n_meshes - 1; ++i)
        {
            pairs_map_distortion.push_back({ i, i + 1 });
        }
    }
    else if (_pairs_mode == DistortionPairs::Cycle)
    {
        for (int i = 0; i < n_meshes; ++i)
        {
            pairs_map_distortion.push_back({ i, (i + 1) % n_meshes });
        }
    }
    else if (_pairs_mode == DistortionPairs::Star)
    {
        for (int i = 1; i < n_meshes; ++i)
        {
            pairs_map_distortion.push_back({ 0, i });
        }
    }
    else
        ISM_ERROR_throw("");

    for (int i = 0; i < (int)pairs_map_distortion.size(); i++)
        ISM_DEBUG_OUT(pairs_map_distortion[i].first << " " << pairs_map_distortion[i].second);
    ISM_DEBUG_VAR(pairs_map_distortion.size());
}

namespace
{

Vec3d
embed_bary(
        const BarycentricPoint& _bary,
        const TriMesh& _mesh,
        const ExternalProperty<VH, Vec3d>& _emb)
{
    return _bary.interpolate(_emb, _mesh).normalized();
}

}

BarycentricPoint map_point(
        const MapState& _map_state,
        const int _idx_from,
        const int _idx_to,
        const BarycentricPoint& _bary_from,
        const TriMesh& _mesh_embedding_T_from,
        const BSPTree& _bsp_T_from)
{
    const TriMesh& mesh_A = _map_state.meshes_input[_idx_from];
    const TriMesh& mesh_B = _map_state.meshes_input[_idx_to];

    // Barycentric point on mesh A --> point in embedding space
    const Vec3d p_emb_A = embed_bary(_bary_from, mesh_A, _map_state.embeddings_input[_idx_from]);

    SFH fh_T;
    double alpha_T, beta_T, gamma_T;
    bsp_tree_barys_face(p_emb_A, _mesh_embedding_T_from, _bsp_T_from, alpha_T, beta_T, gamma_T, fh_T);
    BarycentricPoint bary_T(fh_T, alpha_T, beta_T, _map_state.mesh_T);

    // Barycentric point on mesh T (with embedding T_A) --> Point in embedding space
    const Vec3d p_emb_B = embed_bary(bary_T, _map_state.mesh_T, _map_state.embeddings_T[_idx_to]);

    // Point in embedding space --> Barycentric point on mesh B
    SFH fh_B;
    double alpha_B, beta_B, gamma_B;
    bsp_tree_barys_face(p_emb_B, _map_state.meshes_embeddings_input[_idx_to], _map_state.bsp_embeddings_input[_idx_to], alpha_B, beta_B, gamma_B, fh_B);
    BarycentricPoint bary_B(fh_B, alpha_B, beta_B, mesh_B);
    return bary_B;
}

}
