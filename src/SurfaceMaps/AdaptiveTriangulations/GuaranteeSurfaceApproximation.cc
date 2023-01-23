/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Dörte Pieper, Patrick Schmidt
 */

#include "GuaranteeSurfaceApproximation.hh"

#include <SurfaceMaps/AdaptiveTriangulations/Helpers.hh>
#include <SurfaceMaps/AdaptiveTriangulations/AssignVerticesToFaces.hh>
#include <SurfaceMaps/AdaptiveTriangulations/AdaptiveTargetEdgeLength.hh>
#include <SurfaceMaps/AdaptiveTriangulations/OptimizeMapEnergies.hh>
#include <SurfaceMaps/AdaptiveTriangulations/LiftToSurface.hh>
#include <SurfaceMaps/Misc/ConstantCurvatureGeometry.hh>

namespace SurfaceMaps
{

void split_triangle_on_sphere(
        MapState& _map_state,
        const SFH& _fh_T,
        const std::vector<Vec3d>& _pms)
{
    SVH vh_a, vh_b, vh_c;
    handles(_map_state.mesh_T, _fh_T, vh_a, vh_b, vh_c);

    auto v = _map_state.mesh_T.add_vertex(_pms[0]); // Position does not matter
    ISM_ASSERT_EQ(v.idx(), int(_map_state.mesh_T.n_vertices()) - 1);
    _map_state.mesh_T.split(_fh_T, v);

    for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
    {
        // Set embeddings
        _map_state.embeddings_T[i].resize(_map_state.mesh_T);
        _map_state.embeddings_T[i][v] = _pms[i];

        // Change vertex assignment to new triangles
        std::vector<SVH> before_vhs_i = _map_state.maps_Tf_inputvs[i][_fh_T];
        std::vector<SVH> vhs_mca_i, vhs_mab_i, vhs_mbc_i;

        assign_vertices_to_triangle(before_vhs_i, vhs_mca_i, _map_state.embeddings_T[i][v], _map_state.embeddings_T[i][vh_c], _map_state.embeddings_T[i][vh_a], _map_state.embeddings_input[i]);
        assign_vertices_to_triangle(before_vhs_i, vhs_mab_i, _map_state.embeddings_T[i][v], _map_state.embeddings_T[i][vh_a], _map_state.embeddings_T[i][vh_b], _map_state.embeddings_input[i]);
        assign_vertices_to_triangle(before_vhs_i, vhs_mbc_i, _map_state.embeddings_T[i][v], _map_state.embeddings_T[i][vh_b], _map_state.embeddings_T[i][vh_c], _map_state.embeddings_input[i]);
        // Set
        _map_state.maps_Tf_inputvs[i][_fh_T] = vhs_mca_i;
        _map_state.maps_Tf_inputvs[i].push_back(vhs_mab_i);
        _map_state.maps_Tf_inputvs[i].push_back(vhs_mbc_i);
    }
}

// Distance of point p on mesh A/B to the corresponding point on lifted mesh T
double distance_to_lifted_T(
        const Vec3d& _p_sphere,
        const Vec3d& _p,
        const Vec3d& _a_sphere,
        const Vec3d& _b_sphere,
        const Vec3d& _c_sphere,
        const Vec3d& _a_lifted,
        const Vec3d& _b_lifted,
        const Vec3d& _c_lifted)
{
    // Compute barycentric coordinates (in sphere ambient space) of A/B vertex in T triangle
    double alpha, beta, gamma;
    barys_abc_3d(_a_sphere, _b_sphere, _c_sphere, _p_sphere, Spherical, alpha, beta);
    gamma = 1.0 - alpha - beta;

    // Compute base point of A/B vertex on lifted T triangle
    Vec3d p_base = _a_lifted * alpha + _b_lifted * beta + _c_lifted * gamma;

    return (_p - p_base).norm();
}

double max_dist_verts_T(
        const MapState& _map_state)
{
    double max_dist = 0.0;
    for (auto fh : _map_state.mesh_T.faces())
    {
        VH vh_a, vh_b, vh_c;
        handles(_map_state.mesh_T, fh, vh_a, vh_b, vh_c);

        for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
        {

            // Compute lifted positions of T vertices of face
            const Vec3d a_lifted_i = lift_vertex_to_surface(_map_state.embeddings_T[i][vh_a], _map_state.meshes_input[i], _map_state.meshes_embeddings_input[i], _map_state.bsp_embeddings_input[i]);
            const Vec3d b_lifted_i = lift_vertex_to_surface(_map_state.embeddings_T[i][vh_b], _map_state.meshes_input[i], _map_state.meshes_embeddings_input[i], _map_state.bsp_embeddings_input[i]);
            const Vec3d c_lifted_i = lift_vertex_to_surface(_map_state.embeddings_T[i][vh_c], _map_state.meshes_input[i], _map_state.meshes_embeddings_input[i], _map_state.bsp_embeddings_input[i]);

            // Distances for mesh i for face fh
            for(auto vh : _map_state.maps_Tf_inputvs[i][fh])
            {
                auto dist = distance_to_lifted_T(_map_state.embeddings_input[i][vh], _map_state.meshes_input[i].point(vh),
                                                 _map_state.embeddings_T[i][vh_a], _map_state.embeddings_T[i][vh_b], _map_state.embeddings_T[i][vh_c],
                                                 a_lifted_i, b_lifted_i, c_lifted_i);
                max_dist = fmax(max_dist, dist);

            }
        }
    }

    return max_dist;
}

double n_verts_over_bound(
        const MapState& _map_state,
        const double& _bound)
{
    int n = 0;
    for (auto fh : _map_state.mesh_T.faces())
    {
        VH vh_a, vh_b, vh_c;
        handles(_map_state.mesh_T, fh, vh_a, vh_b, vh_c);

        for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
        {

            // Compute lifted positions of T vertices of face
            const Vec3d a_lifted_i = lift_vertex_to_surface(_map_state.embeddings_T[i][vh_a], _map_state.meshes_input[i], _map_state.meshes_embeddings_input[i], _map_state.bsp_embeddings_input[i]);
            const Vec3d b_lifted_i = lift_vertex_to_surface(_map_state.embeddings_T[i][vh_b], _map_state.meshes_input[i], _map_state.meshes_embeddings_input[i], _map_state.bsp_embeddings_input[i]);
            const Vec3d c_lifted_i = lift_vertex_to_surface(_map_state.embeddings_T[i][vh_c], _map_state.meshes_input[i], _map_state.meshes_embeddings_input[i], _map_state.bsp_embeddings_input[i]);

            // Distances for mesh i for face fh
            for(auto vh : _map_state.maps_Tf_inputvs[i][fh])
            {
                auto dist = distance_to_lifted_T(_map_state.embeddings_input[i][vh], _map_state.meshes_input[i].point(vh),
                                                 _map_state.embeddings_T[i][vh_a], _map_state.embeddings_T[i][vh_b], _map_state.embeddings_T[i][vh_c],
                                                 a_lifted_i, b_lifted_i, c_lifted_i);
                if(dist >= _bound)
                    n++;

            }
        }
    }

    return n;
}

void split_triangle_or_edge_at_p(
        MapState& _map_state,
        const SFH& _fh_T,
        const VH& _max_dist_vertex,
        const int& _max_dist_mesh_idx)
{
    SVH vh_a, vh_b, vh_c;
    handles(_map_state.mesh_T, _fh_T, vh_a, vh_b, vh_c);

    Vec3d a_sphere = _map_state.embeddings_T[_max_dist_mesh_idx][vh_a];
    Vec3d b_sphere = _map_state.embeddings_T[_max_dist_mesh_idx][vh_b];
    Vec3d c_sphere = _map_state.embeddings_T[_max_dist_mesh_idx][vh_c];
    ISM_ASSERT(in_triangle_inclusive_3d(_map_state.embeddings_input[_max_dist_mesh_idx][_max_dist_vertex], a_sphere, b_sphere, c_sphere));

    // Compute barycentric coordinates (in sphere ambient space) of A/B vertex in T triangle
    double alpha, beta, gamma;
    barys_abc_3d(a_sphere, b_sphere, c_sphere, _map_state.embeddings_input[_max_dist_mesh_idx][_max_dist_vertex], Spherical, alpha, beta);
    gamma = 1.0 - alpha - beta;

    // New points on spheres in the triangle
    std::vector<Vec3d> pms; pms.reserve(_map_state.meshes_input.size());
    for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
    {
        pms.push_back(_map_state.embeddings_T[i][vh_a] * alpha + _map_state.embeddings_T[i][vh_b] * beta + _map_state.embeddings_T[i][vh_c] * gamma);
        pms[i].normalize();

        // Test barrier energy, if triangle split would invert -> vertex close to edge -> edge split instead
        ISM_ASSERT_FINITE(barrier_energy(pms[i], _map_state.embeddings_T[i][vh_c], _map_state.embeddings_T[i][vh_a]));
        ISM_ASSERT_FINITE(barrier_energy(pms[i], _map_state.embeddings_T[i][vh_a], _map_state.embeddings_T[i][vh_b]));
        ISM_ASSERT_FINITE(barrier_energy(pms[i], _map_state.embeddings_T[i][vh_b], _map_state.embeddings_T[i][vh_c]));
    }

    // Split triangle
    split_triangle_on_sphere(_map_state, _fh_T, pms);

    //ISM_DEBUG_OUT(_map_state.mesh_T.n_faces());

    // for DEBUG TODO delete
    if(n_flipped_volume(_map_state)>0 && false)
    {
        ISM_WARNING("triangle split leads to flipped triangles");
    }
}


// TODO make sure assignment of map_Tf_Avs and map_Tf_Bvs was done before?
void satisfy_approx_bound_face_splits(
        MapState& _map_state,
        const double& _bound,
        const int &_max_iters)
{
    ISM_ASSERT(sanity_check_assignments(_map_state));
    ISM_ASSERT(sanity_check_all_vertices_assigned(_map_state));

    ExternalProperty<FH, bool> all_verts_under_bound(_map_state.mesh_T, false);

    int cur_iter = 0;
    while(max_dist_verts_T(_map_state) >= _bound && cur_iter < _max_iters)
    {
        cur_iter++;
        for (auto fh : _map_state.mesh_T.faces())
        {
            if(all_verts_under_bound[fh])
                continue;

            VH vh_a, vh_b, vh_c;
            handles(_map_state.mesh_T, fh, vh_a, vh_b, vh_c);

            double max_dist = 0.0;
            int max_dist_mesh_idx;
            VH max_dist_vertex;
            // Find farthest away vertex in the face for all meshes_input
            for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
            {
                // Compute lifted positions of T vertices of face
                const Vec3d a_lifted_i = lift_vertex_to_surface(_map_state.embeddings_T[i][vh_a], _map_state.meshes_input[i], _map_state.meshes_embeddings_input[i], _map_state.bsp_embeddings_input[i]);
                const Vec3d b_lifted_i = lift_vertex_to_surface(_map_state.embeddings_T[i][vh_b], _map_state.meshes_input[i], _map_state.meshes_embeddings_input[i], _map_state.bsp_embeddings_input[i]);
                const Vec3d c_lifted_i = lift_vertex_to_surface(_map_state.embeddings_T[i][vh_c], _map_state.meshes_input[i], _map_state.meshes_embeddings_input[i], _map_state.bsp_embeddings_input[i]);

                // Distances for mesh i for face fh
                for(auto vh : _map_state.maps_Tf_inputvs[i][fh])
                {
                    auto dist = distance_to_lifted_T(_map_state.embeddings_input[i][vh], _map_state.meshes_input[i].point(vh),
                                                     _map_state.embeddings_T[i][vh_a], _map_state.embeddings_T[i][vh_b], _map_state.embeddings_T[i][vh_c],
                                                     a_lifted_i, b_lifted_i, c_lifted_i);
                    if(dist > max_dist)
                    {
                        max_dist = dist;
                        max_dist_mesh_idx = i;
                        max_dist_vertex = vh;
                    }

                }
            }

            // Split at farthest away vertex
            if(max_dist >= _bound)
            {
                split_triangle_or_edge_at_p(_map_state, fh, max_dist_vertex, max_dist_mesh_idx);
                ISM_ASSERT(sanity_check_assignments(_map_state));
                ISM_ASSERT(sanity_check_all_vertices_assigned(_map_state));
            }

            all_verts_under_bound.resize(_map_state.mesh_T, false);
            if(max_dist < _bound)
            {
                all_verts_under_bound[fh] = true;
            }
        }
    }
    max_dist_verts_T(_map_state);
    for (auto fh : _map_state.mesh_T.faces())
    {
        VH vh_a, vh_b, vh_c;
        handles(_map_state.mesh_T, fh, vh_a, vh_b, vh_c);

        for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
        {
            // Compute lifted positions of T vertices of face
            const Vec3d a_lifted_i = lift_vertex_to_surface(_map_state.embeddings_T[i][vh_a], _map_state.meshes_input[i], _map_state.meshes_embeddings_input[i], _map_state.bsp_embeddings_input[i]);
            const Vec3d b_lifted_i = lift_vertex_to_surface(_map_state.embeddings_T[i][vh_b], _map_state.meshes_input[i], _map_state.meshes_embeddings_input[i], _map_state.bsp_embeddings_input[i]);
            const Vec3d c_lifted_i = lift_vertex_to_surface(_map_state.embeddings_T[i][vh_c], _map_state.meshes_input[i], _map_state.meshes_embeddings_input[i], _map_state.bsp_embeddings_input[i]);

            ISM_ASSERT_FINITE(surface_approx_barrier(_map_state.embeddings_T[i][vh_a], _map_state.embeddings_T[i][vh_b], _map_state.embeddings_T[i][vh_c],
                                                    a_lifted_i, b_lifted_i, c_lifted_i, _map_state.maps_Tf_inputvs[i][fh], _map_state.meshes_input[i], _map_state.embeddings_input[i], _map_state.vertex_areas_input[i], _bound));
        }
    }
}

void split_edges(
        MapState& _map_state,
        const std::vector<SEH>& _edges)
{
    // Split edges
    int n_splits = 0;
    for (auto e : _edges)
    {
        // Choose new vertex position at edge midpoint.
        std::vector<Vec3d> pms; pms.reserve(_map_state.meshes_input.size());
        for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
        {
            pms.push_back(0.5 * (_map_state.embeddings_T[i][e.v0()] + _map_state.embeddings_T[i][e.v1()]));
            pms[i].normalize();
        }

        // Split edge
        const Vec3d p_T = 0.5 * (_map_state.mesh_T.point(e.v0()) + _map_state.mesh_T.point(e.v1())); // This position doesn't matter
        const SVH v_new = _map_state.mesh_T.add_vertex(p_T);
        ISM_ASSERT_EQ(v_new.idx(), int(_map_state.mesh_T.n_vertices()) - 1);
        _map_state.mesh_T.split_edge(e, v_new);

        for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
        {
            // Set embeddings
            _map_state.embeddings_T[i].resize(_map_state.mesh_T);
            _map_state.embeddings_T[i][v_new] = pms[i];
        }
        n_splits++;
    }

    ISM_INFO("Split " << n_splits << " edges in local refinement");

    // Re-compute vertex A/B to face T assignments
    assign_vertices_to_T_faces(_map_state);
}

void split_incident_edges(
        MapState& _map_state,
        const ExternalProperty<FH, bool>& _split_fhs)
{
    // Collect edges to split
    std::vector<SEH> split_ehs;
    for (auto e : _map_state.mesh_T.edges())
    {
        if (_split_fhs[e.h0().face()] || _split_fhs[e.h1().face()])
            split_ehs.push_back(e);
    }

    split_edges(_map_state, split_ehs);
}

void split_very_long_edges(
        MapState& _map_state,
        const double _tel_factor = 10.0)
{
    // Collect edges to split
    std::vector<SEH> split_ehs;
    std::vector<TriMesh> lifted_Ts = lifted_meshes_from_mapstate(_map_state);
    for (auto e : _map_state.mesh_T.edges())
    {
        // Pick target edge length at mid point
        double min_tel = INF_DOUBLE;
        for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
        {
            const double tel = target_edge_length(_map_state.embeddings_T[i][e.v0()], _map_state.embeddings_T[i][e.v1()], _map_state.meshes_input[i], _map_state.meshes_embeddings_input[i], _map_state.tels_input[i], _map_state.bsp_embeddings_input[i]);
            min_tel = std::min(min_tel, tel);
        }

        // Find longest instance of this edge in 3D
        double max_el = 0.0;
        for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
            max_el = std::max(max_el, lifted_Ts[i].calc_edge_length(e));

        if (max_el >= _tel_factor * min_tel)
            split_ehs.push_back(e);
    }

    split_edges(_map_state, split_ehs);
}

void local_refine_tel(
        MapState& _map_state,
        const double& _bound,
        const double& _factor,
        const bool _split)
{
    std::vector<ExternalProperty<VH, bool>> refine; refine.reserve(_map_state.meshes_input.size());
    for (int i = 0; i < (int) _map_state.meshes_input.size(); i++)
    {
         refine.push_back(ExternalProperty<VH, bool>(_map_state.meshes_input[i], false));
    }

    ExternalProperty<FH, bool> over_bound(_map_state.mesh_T, false);
    for (auto fh : _map_state.mesh_T.faces())
    {
        VH vh_a, vh_b, vh_c;
        handles(_map_state.mesh_T, fh, vh_a, vh_b, vh_c);

        // Check if any vertex of input meshes mapped into the face fh of T is over bound
        for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
        {

            // Compute lifted positions of T vertices of face
            const Vec3d a_lifted_i = lift_vertex_to_surface(_map_state.embeddings_T[i][vh_a], _map_state.meshes_input[i], _map_state.meshes_embeddings_input[i], _map_state.bsp_embeddings_input[i]);
            const Vec3d b_lifted_i = lift_vertex_to_surface(_map_state.embeddings_T[i][vh_b], _map_state.meshes_input[i], _map_state.meshes_embeddings_input[i], _map_state.bsp_embeddings_input[i]);
            const Vec3d c_lifted_i = lift_vertex_to_surface(_map_state.embeddings_T[i][vh_c], _map_state.meshes_input[i], _map_state.meshes_embeddings_input[i], _map_state.bsp_embeddings_input[i]);

            // Distances for mesh i for face fh
            for(auto vh : _map_state.maps_Tf_inputvs[i][fh])
            {
                auto cur_dist = distance_to_lifted_T(_map_state.embeddings_input[i][vh], _map_state.meshes_input[i].point(vh),
                                                 _map_state.embeddings_T[i][vh_a], _map_state.embeddings_T[i][vh_b], _map_state.embeddings_T[i][vh_c],
                                                 a_lifted_i, b_lifted_i, c_lifted_i);
                if (cur_dist >= _bound)
                {
                    over_bound[fh] = true;
                    break;
                }
            }
        }

        // Refine if over bound
        if (over_bound[fh])
        {
            // Refine tel of input mesh i
            for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
            {
                for(auto vh : _map_state.maps_Tf_inputvs[i][fh])
                {
                    if (!refine[i][vh])
                    {
                        _map_state.tels_input[i][vh] *= _factor;
                        refine[i][vh] = true;
                    }
                }
            }

            // Refine vertices of input mesh in which barycenter of T face lies
            for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
            {
                Vec3d p_center = (_map_state.embeddings_T[i][vh_a] + _map_state.embeddings_T[i][vh_b] + _map_state.embeddings_T[i][vh_c]) * (1.0/3.0);

                SFH f;
                double alpha, beta, gamma;
                bsp_tree_barys_face(p_center, _map_state.meshes_embeddings_input[i], _map_state.bsp_embeddings_input[i], alpha, beta, gamma, f);

                SVH vh_a_i, vh_b_i, vh_c_i;
                handles(_map_state.meshes_input[i], f, vh_a_i, vh_b_i, vh_c_i);

                for (auto vh_A : {vh_a_i, vh_b_i, vh_c_i})
                {
                    if (!refine[i][vh_A])
                    {
                        _map_state.tels_input[i][vh_A] *= _factor;
                        refine[i][vh_A] = true;
                    }
                }

            }

            // Refine tel of vertices of faces in which vertices of T face lie (maximal 9 vertices of each mesh)
            SFH fh_i;
            double alpha, beta, gamma;
            for (auto vh : {vh_a, vh_b, vh_c})
            {
                for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
                {
                    bsp_tree_barys_face(_map_state.embeddings_T[i][vh], _map_state.meshes_embeddings_input[i], _map_state.bsp_embeddings_input[i], alpha, beta, gamma, fh_i);
                    VH vh_a_A, vh_b_A, vh_c_A;
                    handles(_map_state.meshes_input[i], fh_i, vh_a_A, vh_b_A, vh_c_A);
                    for (auto vh_A : {vh_a_A, vh_b_A, vh_c_A})
                    {
                        if (!refine[i][vh_A])
                        {
                            _map_state.tels_input[i][vh_A] *= _factor;
                            refine[i][vh_A] = true;
                        }
                    }
                }
            }
        }
    }

    if (_split)
    {
        split_incident_edges(_map_state, over_bound);
        split_very_long_edges(_map_state);
    }
}

}
