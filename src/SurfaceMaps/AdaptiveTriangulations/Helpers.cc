/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Dörte Pieper, Patrick Schmidt
 */

#include <SurfaceMaps/AdaptiveTriangulations/Helpers.hh>
#include <SurfaceMaps/AdaptiveTriangulations/LiftToSurface.hh>

#include <TinyAD/Utils/Helpers.hh>
#include <SurfaceMaps/Utils/IO.hh>

namespace SurfaceMaps
{

// Return TriMesh of sphere embedding
TriMesh embedding_to_mesh(
        const TriMesh& _mesh_orig,
        const ExternalProperty<VH, Vec3d>& _embedding)
{
    TriMesh embedded_mesh;
    embedded_mesh.assign(_mesh_orig);
    for (auto v : embedded_mesh.vertices())
    {
        embedded_mesh.set_point(v, _embedding[v]);
    }
    return embedded_mesh;
}

std::vector<TriMesh> lifted_meshes_from_mapstate(
        const MapState& _map_state)
{
    const int n_meshes = _map_state.meshes_input.size();

    std::vector<TriMesh> lifted_Ts(n_meshes);
    for (int i = 0; i < n_meshes; ++i)
    {
        lifted_Ts[i] = embedding_to_mesh(
                    _map_state.mesh_T,
                    lift_to_surface(
                        _map_state.meshes_input[i],
                        _map_state.mesh_T,
                        _map_state.embeddings_input[i],
                        _map_state.embeddings_T[i],
                        _map_state.bsp_embeddings_input[i]));
    }

    return lifted_Ts;
}

/// Compute an abitrary tangent vector of the sphere at position p.
Eigen::Vector3d any_orthogonal(
        const Eigen::Vector3d& _p)
{
    // Find coordinate axis spanning the largest angle with _p.
    // Return cross product that of axis with _p
    Eigen::Vector3d tang;
    double min_abs_dot = INFINITY;
    for (const Eigen::Vector3d& ax : { Eigen::Vector3d(1.0, 0.0, 0.0), Eigen::Vector3d(0.0, 1.0, 0.0), Eigen::Vector3d(0.0, 0.0, 1.0) })
    {
        double abs_dot = fabs(_p.dot(ax));
        if (abs_dot < min_abs_dot)
        {
            min_abs_dot = abs_dot;
            tang = ax.cross(_p).normalized();
        }
    }

    return tang;
}

void compute_local_bases(
        const ExternalProperty<VH, Vec3d>& _embedding_input,
        ExternalProperty<VH, Vec3d>& _B1,
        ExternalProperty<VH, Vec3d>& _B2,
        const TriMesh& _mesh_input)
{
    for (auto vh : _mesh_input.vertices())
    {
        _B1[vh] = any_orthogonal(_embedding_input[vh]);
        _B2[vh] = _embedding_input[vh].cross(_B1[vh]);
    }
}


/// Return number of flipped triangles on sphere with volume calculation
int n_flipped_volume(const MapState& _map_state)
{
    int n_flipped = 0;
    for(auto fh : _map_state.mesh_T.faces())
    {
        auto a_idx = fh.halfedge().to();
        auto b_idx = fh.halfedge().next().to();
        auto c_idx = fh.halfedge().from();
        for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
        {
            Eigen::Vector3<double> a_mani = _map_state.embeddings_T[i][a_idx];
            Eigen::Vector3<double> b_mani = _map_state.embeddings_T[i][b_idx];
            Eigen::Vector3<double> c_mani = _map_state.embeddings_T[i][c_idx];

            const double volume = TinyAD::col_mat(a_mani, b_mani, c_mani).determinant();
            if (volume <= 0.0)
            {
                n_flipped++;
                break;
            }
        }
    }
    return n_flipped;
}

double length_diagonal_of_bounding_box(const TriMesh& _mesh)
{
    MatXd V;
    MatXi F;
    mesh_to_matrix(_mesh, V, F);
    // bounding box points
    Eigen::Vector3d min_point = V.colwise().minCoeff();
    Eigen::Vector3d max_point = V.colwise().maxCoeff();

    return (min_point-max_point).norm();
}

}
