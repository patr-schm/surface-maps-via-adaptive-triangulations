/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Dörte Pieper, Patrick Schmidt
 */

#include <SurfaceMaps/AdaptiveTriangulations/LiftToSurface.hh>

#include <SurfaceMaps/AdaptiveTriangulations/Helpers.hh>
#include <SurfaceMaps/Utils/BSPTree.hh>

namespace SurfaceMaps
{

// lifts embedding of _mesh_T to _mesh_imput mesh
ExternalProperty<VH, Vec3d> lift_to_surface(
        const TriMesh& _mesh_input, // A or B
        const TriMesh& _mesh_T,
        const ExternalProperty<VH, Vec3d>& _embedding_input, // embedding_A or embedding_B
        const ExternalProperty<VH, Vec3d>& _embedding_T,
        const BSPTree& _bsp_input)
{
    ExternalProperty<VH, Vec3d> T_on_surface(_mesh_T, Vec3d(0.0, 0.0, 0.0));

    TriMesh mesh_embedding_input = embedding_to_mesh(_mesh_input, _embedding_input);

    // for each vertex of T
    for (auto v:_mesh_T.vertices())
    {
        Vec3d lifted_pos = lift_vertex_to_surface(_embedding_T[v], _mesh_input, mesh_embedding_input, _bsp_input);

        T_on_surface[v] = lifted_pos;
    }

    return T_on_surface;
}

ExternalProperty<VH, Vec3d> lift_to_surface(
        const TriMesh& _mesh_input, // A or B
        const TriMesh& _mesh_T,
        const ExternalProperty<VH, Vec3d>& _embedding_input, // embedding_A or embedding_B
        const ExternalProperty<VH, Vec3d>& _embedding_T)
{
    TriMesh mesh_embedding_input = embedding_to_mesh(_mesh_input, _embedding_input);
    BSPTree bsp_input(mesh_embedding_input);

    return lift_to_surface(_mesh_input, _mesh_T, _embedding_input, _embedding_T, bsp_input);
}

}
