/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt
 */

#include "SphereEmbeddingTet.hh"

#include <TinyAD/Utils/Timer.hh>
#include <SurfaceMaps/Utils/Dijkstra.hh>
#include <SurfaceMaps/Utils/Parametrization.hh>
#include <SurfaceMaps/Misc/ConstantCurvatureGeometry.hh>
#include <SurfaceMaps/Misc/ConstantCurvatureEmbedding.hh>

#include <queue>

namespace SurfaceMaps
{

/// Simple sphere embedding of genus 0 mesh.
/// (1) Pick two neighboring triangles (a, b, c), (c, a, d) and embed them as two sides of a regular tetrahedron.
/// (2) Find a path between the opposite vertices b, d and embed it via constant speed parametrization.
/// (3) Embed all other vertices on the surface of the tetrahedron via harmonic parametrization and project them to the sphere.
bool sphere_embedding_tet(
        const TriMesh& _mesh_L,
        const EH _eh_main,
        ExternalProperty<VH, Vec3d>& _embedding)
{
    ISM_ASSERT_GEQ(_mesh_L.n_vertices(), 4);
    ISM_ASSERT_EQ(genus(_mesh_L), 0);

    _embedding.init(_mesh_L, Vec3d(0.0, 0.0, 0.0));
    ExternalProperty<VH, bool> fixed_vertices(_mesh_L, false);

    // Pick two neighboring triangles and embed vertices as regular tetrahedron
    VH vh_a, vh_b, vh_c, vh_d = VH(-1);

    const HEH heh_to_a = _mesh_L.halfedge_handle(_eh_main, 0);
    const HEH heh_to_b = _mesh_L.next_halfedge_handle(heh_to_a);
    const HEH heh_to_c = _mesh_L.next_halfedge_handle(heh_to_b);
    const HEH heh_to_d = _mesh_L.next_halfedge_handle(_mesh_L.opposite_halfedge_handle(heh_to_a));
    const HEH heh_from_d = _mesh_L.next_halfedge_handle(heh_to_d);

    vh_a = _mesh_L.to_vertex_handle(heh_to_a);
    vh_b = _mesh_L.to_vertex_handle(heh_to_b);
    vh_c = _mesh_L.to_vertex_handle(heh_to_c);
    vh_d = _mesh_L.to_vertex_handle(heh_to_d);

    _embedding[vh_a] = Vec3d(std::sqrt(8.0 / 9.0), 0.0, -1.0 / 3.0);
    _embedding[vh_b] = Vec3d(-std::sqrt(2.0 / 9.0), -std::sqrt(2.0 / 3.0), -1.0 / 3.0);
    _embedding[vh_c] = Vec3d(-std::sqrt(2.0 / 9.0), std::sqrt(2.0 / 3.0), -1.0 / 3.0);
    _embedding[vh_d] = Vec3d(0.0, 0.0, 1.0);
    ISM_ASSERT_EPS(_embedding[vh_a].norm(), 1.0, 1e-9);
    ISM_ASSERT_EPS(_embedding[vh_b].norm(), 1.0, 1e-9);
    ISM_ASSERT_EPS(_embedding[vh_c].norm(), 1.0, 1e-9);
    ISM_ASSERT_EPS(_embedding[vh_d].norm(), 1.0, 1e-9);
    ISM_ASSERT(!flipped_or_degenerate(_mesh_L.face_handle(heh_to_a), _mesh_L, _embedding, Spherical));
    ISM_ASSERT(!flipped_or_degenerate(_mesh_L.opposite_face_handle(heh_to_a), _mesh_L, _embedding, Spherical));
    fixed_vertices[vh_a] = true;
    fixed_vertices[vh_b] = true;
    fixed_vertices[vh_c] = true;
    fixed_vertices[vh_d] = true;

//    {
//        auto v = gv::view();
//        auto style = default_style();
//        view_mesh(_mesh_L);
//    }

//    {
//        auto v = gv::view();
//        auto style = default_style();
//        view_wireframe_on_model(_mesh_L, _embedding, Spherical);
//    }

    // Find path between opposite vertices (b and d)
    ExternalProperty<EH, bool> blocked_edges(_mesh_L, false);

    blocked_edges[_mesh_L.edge_handle(heh_to_a)] = true;
    blocked_edges[_mesh_L.edge_handle(heh_to_b)] = true;
    blocked_edges[_mesh_L.edge_handle(heh_to_c)] = true;
    blocked_edges[_mesh_L.edge_handle(heh_to_d)] = true;
    blocked_edges[_mesh_L.edge_handle(heh_from_d)] = true;

    PrimalPath path = find_primal_path(vh_b, vh_d, _mesh_L, blocked_edges, fixed_vertices);

    // Distribute vertices along the segment (b, d) with unit speed
    {
        ISM_ASSERT(!path.hehs.empty());
        ISM_ASSERT_EQ(_mesh_L.from_vertex_handle(path.hehs.front()), vh_b);
        ISM_ASSERT_EQ(_mesh_L.to_vertex_handle(path.hehs.back()), vh_d);

        double total_length = 0;
        for (const HEH heh : path.hehs)
            total_length += _mesh_L.calc_edge_length(heh);

        double length = 0;
        for (int i = 0; i < (int)path.hehs.size() - 1; ++i)
        {
            const HEH heh = path.hehs[i];
            const VH vh = _mesh_L.to_vertex_handle(heh);

            length += _mesh_L.calc_edge_length(heh);

            const double lambda = length / total_length;
            ISM_ASSERT_G(lambda, 0.0);
            ISM_ASSERT_L(lambda, 1.0);

            _embedding[vh] = (1.0 - lambda) * _embedding[vh_b] + lambda * _embedding[vh_d];
            fixed_vertices[vh] = true;
        }
    }

    // Laplacian smoothing of unfixed vertices
    const int n_iters = 100000;
    int iter = 0;
    for (; iter < n_iters; ++iter)
    {
        const ExternalProperty<VH, Vec3d> prev_embedding = _embedding;

        constexpr double damping = 0.5;
        for (const VH vh : _mesh_L.vertices())
        {
            if (fixed_vertices[vh])
                continue;

            Vec3d avg(0.0, 0.0, 0.0);
            double weight_sum = 0.0;
            for (const HEH heh : _mesh_L.voh_range(vh))
            {
                const VH vh_neigh = _mesh_L.to_vertex_handle(heh);
                const double weight = 1.0;//mean_value_weight(_mesh_L, heh);
                avg += weight * prev_embedding[vh_neigh];
                weight_sum += weight;
            }
            avg /= weight_sum;

            _embedding[vh] += damping * (avg - _embedding[vh]);
        }

        // Converged?
        double step = 0.0;
        for (const VH vh : _mesh_L.vertices())
            step += (_embedding[vh] - prev_embedding[vh]).norm();
        step /= _mesh_L.n_vertices();

        if (step < 1e-12)
            break;
    }

    if (iter < n_iters)
    {
//        ISM_INFO("Harmonic parametrization for sphere embedding converged in " << iter + 1 << " iterations.");
    }
    else
    {
        ISM_WARNING("Harmonic parametrization for sphere embedding did not converge in " << iter << " iterations.");
    }

    // Project to sphere
    for (const VH vh : _mesh_L.vertices())
        _embedding[vh] = _embedding[vh].normalized();

    // Count flips
    int n_flips = 0;
    for (const FH fh : _mesh_L.faces())
    {
        if (flipped_or_degenerate(fh, _mesh_L, _embedding, Spherical))
            ++n_flips;
    }

    if (n_flips > 0)
    {
        ISM_ERROR_throw(n_flips << " flips in layout mesh sphere embedding.");
        return false;
    }

    return true;
}

}
