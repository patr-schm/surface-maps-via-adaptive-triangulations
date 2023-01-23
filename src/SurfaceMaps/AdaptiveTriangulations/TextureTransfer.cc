/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Janis Born, Patrick Schmidt
 */

#include "TextureTransfer.hh"

#include <SurfaceMaps/AdaptiveTriangulations/Helpers.hh>
#include <SurfaceMaps/Utils/Subdivision.hh>
#include <SurfaceMaps/Utils/IO.hh>

namespace SurfaceMaps
{

gv::SharedRenderable make_renderable_smooth(
        const TriMesh& _mesh,
        const glow::SharedTexture2D& _texture)
{
    ISM_ASSERT(_mesh.has_halfedge_texcoords2D());

    for (auto v : _mesh.vertices())
        ISM_ASSERT_FINITE_MAT(_mesh.point(v));

    // Convert to polymesh
    pm::Mesh m;
    auto pos = to_polymesh_tg(_mesh, m);
    auto normals = pm::vertex_normals_uniform(pos);
    auto r = gv::make_renderable(gv::polygons(pos).normals(normals));

    // Set texture coordinates
    auto uvs = m.halfedges().map([&] (auto h)
    {
        const VH vh_from(h.vertex_from().idx.value);
        const VH vh_to(h.vertex_to().idx.value);
        const HEH heh = _mesh.find_halfedge(vh_from, vh_to);
        return tg::pos2(_mesh.texcoord2D(heh)[0], _mesh.texcoord2D(heh)[1]);
    });

    // Flip texture due to different conventions
    configure(*r, gv::textured(uvs, _texture).flip());

    return r;
}

TriMesh transfer_texture_and_subdiv(
        const MapState& _map_state,
        const int _idx_from,
        const int _idx_to,
        const int _n_subdiv)
{
    const TriMesh& mesh_A = _map_state.meshes_input[_idx_from];
    const TriMesh& mesh_B = _map_state.meshes_input[_idx_to];

    // Subdivide mesh A so we can sample the map at an increased resolution
    TriMesh mesh_A_subdiv = mesh_A;
    RefinementMap ref_map = identity_refinement_map(mesh_A_subdiv);
    for (int subdiv_i = 0; subdiv_i < _n_subdiv; ++subdiv_i)
        ref_map = subdivide_1_to_4(mesh_A_subdiv, ref_map);

    // Map refined positions of A onto B using the computed map
    TriMesh mesh_A_subdiv_on_B = mesh_A_subdiv;
    TriMesh mesh_embedding_T_A = embedding_to_mesh(_map_state.mesh_T, _map_state.embeddings_T[_idx_from]);
    BSPTree bsp_T_X = BSPTree(mesh_embedding_T_A);
    for (const auto& vh : mesh_A_subdiv_on_B.vertices())
    {
        const HEH& in_heh = vh.in(); // Positions are continuous at vertices, so we can pick any incoming halfedge.
        const BarycentricPoint bary_A = ref_map[in_heh];
        const BarycentricPoint bary_B = map_point(_map_state, _idx_from, _idx_to, bary_A, mesh_embedding_T_A, bsp_T_X);
        mesh_A_subdiv_on_B.point(vh) = bary_B.point(mesh_B);
    }

    // Transfer texture coordinates from coarse to subdivided mesh
    refine_tex_coords(mesh_A, mesh_A_subdiv_on_B, ref_map);

    return mesh_A_subdiv_on_B;
}

}
