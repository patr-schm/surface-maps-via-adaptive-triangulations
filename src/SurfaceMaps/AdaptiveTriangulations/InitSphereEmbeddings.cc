/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt
 */

#include "InitSphereEmbeddings.hh"

#include <TinyAD/Utils/Timer.hh>
#include <SurfaceMaps/Utils/Genus.hh>
#include <SurfaceMaps/Utils/IO.hh>
#include <SurfaceMaps/Utils/MeshNormalization.hh>

namespace SurfaceMaps
{

bool init_map(
        MapState& _map_state,
        const std::vector<fs::path>& _mesh_paths,
        const std::vector<fs::path>& _landmark_paths,
        const std::vector<fs::path>& _embedding_paths,
        const bool _align_meshes)
{
    ISM_ASSERT_EQ(_mesh_paths.size(), _landmark_paths.size());
    ISM_ASSERT_EQ(_mesh_paths.size(), _embedding_paths.size());
    ISM_ASSERT_EQ(_map_state.meshes_input.size(), 0);

    // Load meshes, load landmarks, normalize meshes
    for (int i = 0; i < (int)_mesh_paths.size(); ++i)
    {
        _map_state.meshes_input.push_back(read_mesh(_mesh_paths[i]));
        _map_state.landmarks_input.push_back(read_landmarks(_landmark_paths[i]));
        ISM_ASSERT_EQ(_map_state.landmarks_input[i].size(), _map_state.landmarks_input[0].size());

        if (!closed_surface(_map_state.meshes_input[i]) || genus(_map_state.meshes_input[i]) != 0)
        {
            ISM_WARNING("Input mesh not genus 0");
            return false;
        }

        center_mesh(_map_state.meshes_input[i]);
        normalize_surface_area(_map_state.meshes_input[i]);

        // Align all meshes to first
        if (i != 0 && _align_meshes)
        {
            if (_map_state.landmarks_input[i].size() < 4)
                ISM_DEBUG_OUT("Cannot align meshes: Too few landmarks.")
            else
                align_rigid(_map_state.landmarks_input[0], _map_state.landmarks_input[i], _map_state.meshes_input[0], _map_state.meshes_input[i]);
        }
    }

    // Write normalized input meshes
    for (int i = 0; i < (int)_mesh_paths.size(); ++i)
    {
        const fs::path mesh_path = _embedding_paths[i].parent_path() / ("input_normalized_" + pad_integer(i, 2) + ".obj");
        write_mesh(_map_state.meshes_input[i], mesh_path);
    }

    // Load or compute sphere embeddings
    for (int i = 0; i < (int)_mesh_paths.size(); ++i)
    {
        if (!fs::exists(_embedding_paths[i]))
        {
            ExternalProperty<VH, Vec3d> embedding_tmp = multi_res_sphere_embedding(_map_state.meshes_input[i]);
            ISM_ASSERT(sphere_embedding_bijective(_map_state.meshes_input[i], embedding_tmp));

            // Rotation-align all sphere embeddings to first
            if (i != 0)
            {
                align_rotation(
                            _map_state.landmarks_input[0], _map_state.landmarks_input[i],
                            _map_state.meshes_input[0], _map_state.meshes_input[i],
                            _map_state.embeddings_input[0], embedding_tmp);
            }

            // Save sphere embedding
            write_embedding(_map_state.meshes_input[i], embedding_tmp, _embedding_paths[i]);
        }

        // Always load sphere embedding from file (to have same numerics)
        _map_state.embeddings_input.push_back(read_embedding(_map_state.meshes_input[i], _embedding_paths[i]));
    }

    // Init map
    _map_state.compute_bsp_trees();

    // Init common triangulation with mesh_A
    _map_state.mesh_T = _map_state.meshes_input[0];
    _map_state.landmarks_T = _map_state.landmarks_input[0];
    for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
        _map_state.embeddings_T.push_back(_map_state.embeddings_input[0]);

    _map_state.set_distortion_pairs(DistortionPairs::Star);

    return true;
}

void compute_input_sphere_embeddings(
        const TriMesh& _mesh_A,
        const TriMesh& _mesh_B,
        const std::vector<VH>& _landmarks_A,
        const std::vector<VH>& _landmarks_B,
        ExternalProperty<VH, Vec3d>& _embedding_A,
        ExternalProperty<VH, Vec3d>& _embedding_B,
        const MultiResSphereEmbeddingSettings& _settings)
{
    TinyAD::Timer timer(__FUNCTION__);

    // Compute bijective sphere embeddings
    _embedding_A = multi_res_sphere_embedding(_mesh_A, _settings);
    _embedding_B = multi_res_sphere_embedding(_mesh_B, _settings);
    ISM_ASSERT(sphere_embedding_bijective(_mesh_A, _embedding_A));
    ISM_ASSERT(sphere_embedding_bijective(_mesh_B, _embedding_B));

    // Compute optimal rotation w.r.t. landmarks on sphere. (Rotates B)
    align_rotation(_landmarks_A, _landmarks_B, _mesh_A, _mesh_B, _embedding_A, _embedding_B);
}

}
