/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Dörte Pieper, Patrick Schmidt
 */

#include <SurfaceMaps/AdaptiveTriangulations/Remeshing.hh>

#include <SurfaceMaps/AdaptiveTriangulations/Helpers.hh>
#include <SurfaceMaps/AdaptiveTriangulations/OptimizeVertex.hh>
#include <SurfaceMaps/AdaptiveTriangulations/OptimizeMapEnergies.hh>
#include <SurfaceMaps/AdaptiveTriangulations/AssignVerticesToFaces.hh>
#include <SurfaceMaps/Misc/ConstantCurvatureGeometry.hh>
#include <SurfaceMaps/Utils/GarbageCollection.hh>
#include <SurfaceMaps/Utils/Helpers.hh>
#include <TinyAD/Utils/Timer.hh>
#include <SurfaceMaps/Utils/IO.hh>
#include <TinyAD/Utils/Helpers.hh>
#include <TinyAD/Utils/Timer.hh>

namespace SurfaceMaps
{

namespace
{

/// Return if triangles after edge split are injective
bool injective_after_split(
        const MapState& _map_state,
        const SEH& _eh,
        const std::vector<Vec3d>& p_embs)
{
    const auto v_from = _eh.h0().from();
    const auto v_to = _eh.h0().to();
    const auto v_1 = _eh.h0().next().to();
    const auto v_2 = _eh.h0().opp().next().to();

    for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
    {
        if (!ccw_exclusive_3d(_map_state.embeddings_T[i][v_from], p_embs[i], _map_state.embeddings_T[i][v_1]))
            return false;
        if (!ccw_exclusive_3d(_map_state.embeddings_T[i][v_1], p_embs[i], _map_state.embeddings_T[i][v_to]))
            return false;
        if (!ccw_exclusive_3d(_map_state.embeddings_T[i][v_2], p_embs[i], _map_state.embeddings_T[i][v_from]))
            return false;
        if (!ccw_exclusive_3d(_map_state.embeddings_T[i][v_to], p_embs[i], _map_state.embeddings_T[i][v_2]))
            return false;
    }
    return true;
}

/// Compute the improvement in energy that the split of the edge _eh would yield
double split_energy_improvement(
        const MapState& _map_state,
        const SEH& _eh,
        const AdaptiveTriangulationsSettings& _settings,
        std::vector<std::vector<Triangle>>& _tris_after) // needed outside to update assignments in case of split
{
    /*        +  d
     *       /|\
     *      / | \
     *     /  |m \
     *  a +---+---+ c
     *     \  |  /
     *      \ | /
     *       \|/
     *        +  b
     */
    const SVH va = _eh.h0().to();
    const SVH vb = _eh.h0().next().to();
    const SVH vc = _eh.h0().from();
    const SVH vd = _eh.h0().opp().next().to();
    const SFH f_l = _eh.h0().face();       // before split abc, after split mab
    const SFH f_r = _eh.h0().opp().face(); // before split acd, after split mda
    const SFH f_l_new(_map_state.mesh_T.n_faces(), &_map_state.mesh_T); // after split mbc
    const SFH f_r_new(_map_state.mesh_T.n_faces() + 1, &_map_state.mesh_T); // after split mcd

    // Choose new vertex position at edge midpoint.
    std::vector<Vec3d> pms; pms.reserve(_map_state.meshes_input.size());
    for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
    {
        pms.push_back(0.5 * (_map_state.embeddings_T[i][_eh.v0()] + _map_state.embeddings_T[i][_eh.v1()]));
        pms[i].normalize();
    }

    // Triangles before split
    std::vector<std::vector<Triangle>> tris_before; tris_before.reserve(_map_state.meshes_input.size());
    for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
    {
        tris_before.push_back({
                                  Triangle(f_l, _map_state.mesh_T, _map_state.embeddings_T[i], _map_state.maps_Tf_inputvs[i], _settings), // abc_i
                                  Triangle(f_r, _map_state.mesh_T, _map_state.embeddings_T[i], _map_state.maps_Tf_inputvs[i], _settings), // acd_i
                              });
    }

    // Triangles after split (variable vertex always has to be first)
    _tris_after.reserve(_map_state.meshes_input.size());
    for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
    {
        _tris_after.push_back({
                                 Triangle(f_l, pms[i], _map_state.embeddings_T[i][va], _map_state.embeddings_T[i][vb]), // mab_A
                                 Triangle(f_l_new, pms[i], _map_state.embeddings_T[i][vb], _map_state.embeddings_T[i][vc]), // mbc_A
                                 Triangle(f_r_new, pms[i], _map_state.embeddings_T[i][vc], _map_state.embeddings_T[i][vd]), // mcd_A
                                 Triangle(f_r, pms[i], _map_state.embeddings_T[i][vd], _map_state.embeddings_T[i][va]), // mda_A
                             });
    }
    ISM_ASSERT_EQ(tris_before.size(), _map_state.meshes_input.size());
    ISM_ASSERT_EQ(_tris_after.size(), _map_state.meshes_input.size());

    // Compute energy before and after. Assign input meshes vertices to after state.
    std::vector<std::vector<SVH>> changing_vhs; changing_vhs.reserve(_map_state.meshes_input.size());
    for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
    {
        changing_vhs.push_back(collect_assigned_vertices(tris_before[i], _settings));
    }

    const double energy_after = eval_and_assign_triangles(_tris_after, changing_vhs, _map_state, _settings);
    ISM_ASSERT_NOT_NAN(energy_after);

    if (energy_after == INFINITY)
        return -INFINITY;
    const double energy_before = eval_triangles(tris_before, _map_state, _settings);

    return energy_before - energy_after;
}

/// Return true if split happened
int split_if_improves_energy(
        MapState& _map_state,
        const AdaptiveTriangulationsSettings& _settings)
{
    /*        +  d
     *       /|\
     *      / | \
     *     /  |m \
     *  a +---+---+ c
     *     \  |  /
     *      \ | /
     *       \|/
     *        +  b
     */
    std::vector<double> energy_improvement(_map_state.mesh_T.n_edges());
    // compute energy improvement in parallel for all edges
    #pragma omp parallel for
    for (int idx_edge = 0; idx_edge < (int)_map_state.mesh_T.n_edges(); idx_edge++)
    {
        SEH eh = SEH(idx_edge, &_map_state.mesh_T);
        std::vector<std::vector<Triangle>> tris_after;
        energy_improvement[idx_edge] = split_energy_improvement(_map_state, eh, _settings, tris_after);
    }

    // Collect all edges with energy improvement.
    std::vector<std::pair<double, int>> improving_edges;
    for (int idx_edge = 0; idx_edge < (int)_map_state.mesh_T.n_edges(); idx_edge++)
    {
        if (energy_improvement[idx_edge] > 0.0)
            improving_edges.push_back({energy_improvement[idx_edge], idx_edge});
    }

    // Sort by improvement in descending order to flip edge with most improvements first
    std::sort(improving_edges.rbegin(), improving_edges.rend());

    int n_splits = 0;
    // Check edges if flip still improves energy and flip if that is the case
    for (auto pair_edge : improving_edges)
    {
        SEH eh = SEH(pair_edge.second, &_map_state.mesh_T);
        //ISM_DEBUG_VAR(pair_edge.first << " " << eh);
        std::vector<std::vector<Triangle>> tris_after;
        double cur_improvement = split_energy_improvement(_map_state, eh, _settings, tris_after);

        // Split if improves energy
        if (cur_improvement > 0.0)
        {
            const SFH f_l_new(_map_state.mesh_T.n_faces(), &_map_state.mesh_T); // after split mbc
            const SFH f_r_new(_map_state.mesh_T.n_faces() + 1, &_map_state.mesh_T); // after split mcd

            // Choose new vertex position at edge midpoint.
            std::vector<Vec3d> pms; pms.reserve(_map_state.meshes_input.size());
            for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
            {
                pms.push_back(0.5 * (_map_state.embeddings_T[i][eh.v0()] + _map_state.embeddings_T[i][eh.v1()]));
                pms[i].normalize();
            }

            // Split edge
            ISM_ASSERT(injective_after_split(_map_state, eh, pms)); // Should not happen because barrier energy would be inf
            const Vec3d p_T = 0.5 * (_map_state.mesh_T.point(eh.v0()) + _map_state.mesh_T.point(eh.v1())); // This position doesn't matter
            const SVH v_new = _map_state.mesh_T.add_vertex(p_T);
            ISM_ASSERT_EQ(v_new.idx(), int(_map_state.mesh_T.n_vertices()) - 1);
            _map_state.mesh_T.split_edge(eh, v_new);

            // Check new face indices
            ISM_ASSERT_EQ(eh.h0().prev().opp().face(), f_l_new);
            ISM_ASSERT_EQ(eh.h1().next().opp().face(), f_r_new);


            for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
            {
                // Set embeddings
                _map_state.embeddings_T[i].resize(_map_state.mesh_T);
                _map_state.embeddings_T[i][v_new] = pms[i];

                // Set new vertex assignments
                _map_state.maps_Tf_inputvs[i].resize(_map_state.mesh_T);
                for (int j = 0; j < (int)tris_after[i].size(); ++j)
                {
                    _map_state.maps_Tf_inputvs[i][tris_after[i][j].fh] = tris_after[i][j].assigned_vhs;
                }
            }
            n_splits++;
        }
    }
    return n_splits;
}


/// Return true if triangles after collapse are injective
bool injective_after_collapse(
        const MapState& _map_state,
        const SHEH& _heh,
        std::vector<Vec3d>& p_embs)
{
    const auto v_from = _heh.from();
    const auto v_to = _heh.to();
    auto cur_heh = _heh.next().opp();

    // Faces around v_to that will change
    while (cur_heh.next().opp().idx() != _heh.idx())
    {
        auto v_a = cur_heh.from();
        ISM_ASSERT_EQ(cur_heh.to().idx(), v_to.idx());
        cur_heh = cur_heh.next();
        auto v_c = cur_heh.to();
        ISM_ASSERT_EQ(cur_heh.from().idx(), v_to.idx());
        cur_heh = cur_heh.opp();
        // check if after collapse counterclockwise
        for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
        {
            if (!ccw_exclusive_3d(_map_state.embeddings_T[i][v_a], p_embs[i], _map_state.embeddings_T[i][v_c]))
                return false;
        }
    }
    cur_heh = cur_heh.prev().opp();

    // Faces around v_from that will change
    while (cur_heh.next().idx() != _heh.idx())
    {
        auto v_a = cur_heh.from();
        ISM_ASSERT_EQ(cur_heh.to().idx(), v_from.idx());
        cur_heh = cur_heh.next();
        auto v_c = cur_heh.to();
        ISM_ASSERT_EQ(cur_heh.from().idx(), v_from.idx());
        cur_heh = cur_heh.opp();
        // check if after collapse counterclockwise
        for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
        {
            if (!ccw_exclusive_3d(_map_state.embeddings_T[i][v_a], p_embs[i], _map_state.embeddings_T[i][v_c]))
                return false;
        }
    }

    return true;
}

/// Compute the improvement in energy that the collapse of the edge _eh would yield
double collapse_energy_improvement(
        const MapState& _map_state,
        const SEH& _eh,
        const AdaptiveTriangulationsSettings& _settings,
        std::vector<std::vector<Triangle>>& _tris_after) // needed outside to update assignments in case of collapse
{
    const auto heh_collapse = _eh.h0();
    const auto v_from = heh_collapse.from();
    const auto v_to = heh_collapse.to();

    // Move remaining vertex to collapsed edge midpoint
    std::vector<Vec3d> pms; pms.reserve(_map_state.meshes_input.size());
    for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
    {
        pms.push_back(0.5 * (_map_state.embeddings_T[i][v_from] + _map_state.embeddings_T[i][_eh.v1()]));
        pms[i].normalize();
    }

    // Skip edge if already deleted
    if (_eh.deleted())
        return -INFINITY;

    // Does connectivity allow a collapse?
    // Skip this check here because OpenMesh can't do this in parallel.
    // Make sure to check before actually performing the collapse.
//    if (!_map_state.mesh_T.is_collapse_ok(heh_collapse))
//        return -INFINITY;

    // Don't collapse edge if incident to landmark vertex
    if (find(_map_state.landmarks_T.begin(), _map_state.landmarks_T.end(), v_from) != _map_state.landmarks_T.end() ||
        find(_map_state.landmarks_T.begin(), _map_state.landmarks_T.end(), v_to) != _map_state.landmarks_T.end())
    {
        return -INFINITY;
    }

    // Faces that will be deleted
    const SFH f_collapse_l = heh_collapse.face();
    const SFH f_collapse_r = heh_collapse.opp().face();

    // Triangles before split
    std::vector<std::vector<Triangle>> tris_before; tris_before.reserve(_map_state.meshes_input.size());
    for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
    {
        std::vector<Triangle> tris_before_i; tris_before_i.reserve(10);
        for (auto h : v_to.outgoing_halfedges())
        {
            tris_before_i.push_back(Triangle(h.face(), _map_state.mesh_T, _map_state.embeddings_T[i], _map_state.maps_Tf_inputvs[i], _settings));
        }
        for (auto h : v_from.outgoing_halfedges())
        {
            SFH f = h.face();
            if (f == f_collapse_l || f == f_collapse_r) // not those faces again
                continue;

            tris_before_i.push_back(Triangle(h.face(), _map_state.mesh_T, _map_state.embeddings_T[i], _map_state.maps_Tf_inputvs[i], _settings));
        }
        tris_before.push_back(tris_before_i);
    }

    _tris_after.reserve(_map_state.meshes_input.size());
    for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
    {
        // Triangles after split
        std::vector<Triangle> tris_after_i; tris_after_i.reserve(8);
        for (auto h : v_to.outgoing_halfedges())
        {
            // Skip deleted faces
            if (h.face() == f_collapse_l || h.face() == f_collapse_r)
                continue;

            // Center vertex has to be the first one!
            tris_after_i.push_back(Triangle(h.face(), pms[i], _map_state.embeddings_T[i][h.to()], _map_state.embeddings_T[i][h.next().to()]));
        }
        for (auto h : v_from.outgoing_halfedges())
        {
            // Skip deleted faces
            if (h.face() == f_collapse_l || h.face() == f_collapse_r)
                continue;

            // Center vertex has to be the first one!
            tris_after_i.push_back(Triangle(h.face(), pms[i], _map_state.embeddings_T[i][h.to()], _map_state.embeddings_T[i][h.next().to()]));
        }
        _tris_after.push_back(tris_after_i);
    }

    // Compute energy before and after. Assign A/B vertices to after state.
    std::vector<std::vector<SVH>> changing_vhs; changing_vhs.reserve(_map_state.meshes_input.size());
    for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
    {
        changing_vhs.push_back(collect_assigned_vertices(tris_before[i], _settings));
    }

    const double energy_after = eval_and_assign_triangles(_tris_after, changing_vhs, _map_state, _settings);
    if (energy_after == INFINITY)
        return -INFINITY;
    const double energy_before = eval_triangles(tris_before, _map_state, _settings);

    return energy_before - energy_after;
}

/// Returns number of collapses
int collapse_if_improves_energy(
        MapState& _map_state,
        const AdaptiveTriangulationsSettings& _settings)
{
    std::vector<double> energy_improvement(_map_state.mesh_T.n_edges());

    // compute energy improvement in parallel for all edges
    #pragma omp parallel for schedule(dynamic)
    for (int idx_edge = 0; idx_edge < (int)_map_state.mesh_T.n_edges(); idx_edge++)
    {
        SEH eh = SEH(idx_edge, &_map_state.mesh_T);
        std::vector<std::vector<Triangle>> tris_after;
        energy_improvement[idx_edge] = collapse_energy_improvement(_map_state, eh, _settings, tris_after);
    }

    // Collect all edges with energy improvement.
    std::vector<std::pair<double, int>> improving_edges;
    for (int idx_edge = 0; idx_edge < (int)_map_state.mesh_T.n_edges(); idx_edge++)
    {
        if (energy_improvement[idx_edge] > 0.0)
            improving_edges.push_back({energy_improvement[idx_edge], idx_edge});
    }

    // Sort by improvement in descending order to flip edge with most improvements first
    std::sort(improving_edges.rbegin(), improving_edges.rend());

    int n_collapses = 0;
    // Check edges if collapse still improves energy and collapse if that is the case
    for (auto pair_edge : improving_edges)
    {
        SEH eh = SEH(pair_edge.second, &_map_state.mesh_T);

        // Does connectivity allow this collapse?
        // This check is skipped in collapse_energy_improvement
        // because OpenMesh cannot check this in parallel.
        if (!_map_state.mesh_T.is_collapse_ok(eh.h0()))
            continue;

        //ISM_DEBUG_VAR(pair_edge.first << " " << eh);
        std::vector<std::vector<Triangle>> tris_after;
        double cur_improvement = collapse_energy_improvement(_map_state, eh, _settings, tris_after);

        // Collapse if improves energy
        if (cur_improvement > 0.0)
        {
            const auto heh_collapse = eh.h0();
            const auto v_from = heh_collapse.from();
            const auto v_to = heh_collapse.to();
            // Move remaining vertex to collapsed edge midpoint
            std::vector<Vec3d> pms; pms.reserve(_map_state.meshes_input.size());
            for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
            {
                pms.push_back(0.5 * (_map_state.embeddings_T[i][v_from] + _map_state.embeddings_T[i][eh.v1()]));
                pms[i].normalize();
            }

            // Collapse edge
            const Vec3d p_T = 0.5 * (_map_state.mesh_T.point(v_from) + _map_state.mesh_T.point(v_to)); // This position does not matter
            ISM_ASSERT(injective_after_collapse(_map_state, heh_collapse, pms)) // Flip should not happen because barrier energy inf
            _map_state.mesh_T.collapse(heh_collapse);
            _map_state.mesh_T.point(v_to) = p_T;

            for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
            {
                // Relocate remaining vertex
                _map_state.embeddings_T[i][v_to] = pms[i];

                // Set new vertex assignments
                for (int j = 0; j < (int)tris_after[i].size(); ++j)
                {
                    _map_state.maps_Tf_inputvs[i][tris_after[i][j].fh] = tris_after[i][j].assigned_vhs;
                }

                // Faces that will be deleted
                const SFH f_collapse_l = heh_collapse.face();
                // Remove vertex assignments for deleted faces (will be removed in garbage collection)
                _map_state.maps_Tf_inputvs[i][f_collapse_l] = {};
            }
            n_collapses++;
        }

    }

    return n_collapses;
}

/// Return true if triangles after edge flip are injective
bool injective_after_flip(
        const MapState& _map_state,
        const SEH& _eh)
{
    const auto heh = _eh.h0();
    const auto v_from = heh.from();
    const auto v_to = heh.to();
    const auto v_next = heh.next().to();
    const auto v_opp = heh.opp().next().to();

    for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
    {
        if (!ccw_exclusive_3d(_map_state.embeddings_T[i][v_to], _map_state.embeddings_T[i][v_next], _map_state.embeddings_T[i][v_opp]))
            return false;
        if (!ccw_exclusive_3d(_map_state.embeddings_T[i][v_from], _map_state.embeddings_T[i][v_opp], _map_state.embeddings_T[i][v_next]))
            return false;
    }

    return true;
}

/// Return if volume spanned by vertices for given embedding on sphere and sphere center is negative (or smaller some value)
bool face_is_flipped_volume(
        const ExternalProperty<VH, Vec3d>& embedding,
        const VH& v_1,
        const VH& v_2,
        const VH& v_3)
{
    return TinyAD::col_mat(embedding[v_1], embedding[v_2], embedding[v_3]).determinant() <=  1e-16;
}

/// Return if volume spanned by triangles after flip and spherecenter would be positive -> not flipped on sphere
bool is_flip_ok_volume(
        const MapState& _map_state,
        const SEH& _eh)
{
    const auto heh = _eh.h0();
    const auto v_from = heh.from();
    const auto v_to = heh.to();
    const auto v_next = heh.next().to();
    const auto v_opp = heh.opp().next().to();

    for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
    {
        if (face_is_flipped_volume(_map_state.embeddings_T[i], v_to, v_next, v_opp))
            return false;
        if (face_is_flipped_volume(_map_state.embeddings_T[i], v_from, v_opp, v_next))
            return false;
    }


    return true;
}

/// Compute the improvement in energy that the flip of the edge _eh would yield
double flip_energy_improvement(
        const MapState& _map_state,
        const SEH& _eh,
        const AdaptiveTriangulationsSettings& _settings,
        std::vector<std::vector<Triangle>>& _tris_after) // needed outside to update assignments in case of flip
{
    /*        +  d
     *       / \
     *      /   \
     *     /  f2 \
     *  a +<------+ c
     *     \  f1 /
     *      \   /
     *       \ /
     *        +  b
     */

    const SHEH heh = _eh.h0();
    const SVH va = heh.to();
    const SVH vb = heh.next().to();
    const SVH vc = heh.from();
    const SVH vd = heh.opp().next().to();
    SFH f1 = heh.face();       // before flip: abc, after flip bcd
    SFH f2 = heh.opp().face(); // before flip: acd, after flip abd

    // Check connectivity
    if (!_map_state.mesh_T.is_flip_ok(_eh))
       return -INFINITY;

    // Do not directly connect two landmarks
    if (find(_map_state.landmarks_T.begin(), _map_state.landmarks_T.end(), vb) != _map_state.landmarks_T.end() &&
       find(_map_state.landmarks_T.begin(), _map_state.landmarks_T.end(), vd) != _map_state.landmarks_T.end())
       return -INFINITY;

    // Triangles before and after edge flip
    std::vector<std::vector<Triangle>> tris_before; tris_before.reserve(_map_state.meshes_input.size());
    _tris_after.reserve(_map_state.meshes_input.size());
    for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
    {
        tris_before.push_back({
                                  Triangle(f1, _map_state.mesh_T, _map_state.embeddings_T[i], _map_state.maps_Tf_inputvs[i], _settings), // abc_i
                                  Triangle(f2, _map_state.mesh_T, _map_state.embeddings_T[i], _map_state.maps_Tf_inputvs[i], _settings), // acd_i
                              });
    }
    for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
    {
        _tris_after.push_back({
                                 Triangle(f1, _map_state.embeddings_T[i][vb], _map_state.embeddings_T[i][vc], _map_state.embeddings_T[i][vd]), // bcd_i
                                 Triangle(f2, _map_state.embeddings_T[i][va], _map_state.embeddings_T[i][vb], _map_state.embeddings_T[i][vd]), // abd_A
                             });
    }
    ISM_ASSERT_EQ(tris_before.size(), _map_state.meshes_input.size());
    ISM_ASSERT_EQ(_tris_after.size(), _map_state.meshes_input.size());


    // Compute energy before and after. Assign A/B vertices to after state.
    std::vector<std::vector<SVH>> changing_vhs; changing_vhs.reserve(_map_state.meshes_input.size());
    for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
    {
        changing_vhs.push_back(collect_assigned_vertices(tris_before[i], _settings));
    }

    const double energy_after = eval_and_assign_triangles(_tris_after, changing_vhs, _map_state, _settings);
    if (energy_after == INFINITY)
        return -INFINITY;
    const double energy_before = eval_triangles(tris_before, _map_state, _settings);

    return energy_before - energy_after;
}

/// Returns number of flips
int flip_if_improves_energy(
        MapState& _map_state,
        //const SEH _eh,
        const AdaptiveTriangulationsSettings& _settings)
{
    /*        +  d
     *       / \
     *      /   \
     *     /  f2 \
     *  a +<------+ c
     *     \  f1 /
     *      \   /
     *       \ /
     *        +  b
     */

    std::vector<double> energy_improvement(_map_state.mesh_T.n_edges());
    // compute energy improvement in parallel for all edges
    #pragma omp parallel for
    for (int idx_edge = 0; idx_edge < (int)_map_state.mesh_T.n_edges(); idx_edge++)
    {
        SEH eh = SEH(idx_edge, &_map_state.mesh_T);
        std::vector<std::vector<Triangle>> tris_after;
        energy_improvement[idx_edge] = flip_energy_improvement(_map_state, eh, _settings, tris_after);
    }

    // Collect all edges with energy improvement.
    std::vector<std::pair<double, int>> improving_edges;
    for (int idx_edge = 0; idx_edge < (int)_map_state.mesh_T.n_edges(); idx_edge++)
    {
        if (energy_improvement[idx_edge] > 0.0)
            improving_edges.push_back({energy_improvement[idx_edge], idx_edge});
    }

    // Sort by improvement in descending order to flip edge with most improvements first
    std::sort(improving_edges.rbegin(), improving_edges.rend());

    int n_flips = 0;
    // Check edges if flip still improves energy and flip if that is the case
    for (auto pair_edge : improving_edges)
    {
        SEH eh = SEH(pair_edge.second, &_map_state.mesh_T);
        //ISM_DEBUG_VAR(pair_edge.first << " " << eh);
        std::vector<std::vector<Triangle>> tris_after;
        double cur_improvement = flip_energy_improvement(_map_state, eh, _settings, tris_after);

        // Flip if improves energy
        if (cur_improvement > 0.0)
        {
            if (!is_flip_ok_volume(_map_state, eh)) // Check with some epsilon. Seems to be necessary.
                continue;
                //return false;

            ISM_ASSERT(injective_after_flip(_map_state, eh));

            // Flip
            _map_state.mesh_T.flip(eh);

            // Set new vertex assignments
            for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
            {
                for (int j = 0; j < (int)tris_after[i].size(); ++j)
                {
                    _map_state.maps_Tf_inputvs[i][tris_after[i][j].fh] = tris_after[i][j].assigned_vhs;
                }
            }

            n_flips++;
        }
    }
    return n_flips;

}

}

bool remesh_T(
        MapState& _map_state,
        const AdaptiveTriangulationsSettings& _settings,
        const fs::path& _csv_path,
        std::function<void(const std::string&)> _callback)
{
    // Counters
    int n_splits = 0;
    int n_collapses = 0;
    int n_flips = 0;

    ISM_ASSERT(_settings.w_approx <= 0.0 || sanity_check_all_vertices_assigned(_map_state));
    ISM_ASSERT(_settings.w_approx <= 0.0 || sanity_check_assignments(_map_state));

    // Split if improves energy
    TinyAD::Timer timer_split("timer_split");
    if (_settings.allow_splits)
        n_splits = split_if_improves_energy(_map_state, _settings);
    timer_split.stop();

    ISM_DEBUG_OUT("splits done");

    // Callback
    _callback("after_splits");

    ISM_ASSERT(_settings.w_approx <= 0.0 || sanity_check_all_vertices_assigned(_map_state));
    ISM_ASSERT(_settings.w_approx <= 0.0 || sanity_check_assignments(_map_state));

    // Request status
    _map_state.mesh_T.request_vertex_status();
    _map_state.mesh_T.request_halfedge_status();
    _map_state.mesh_T.request_edge_status();
    _map_state.mesh_T.request_face_status();

    // Collapse edges if improves energy
    TinyAD::Timer timer_collapse("timer collapse");
    if (_settings.allow_collapses)
        n_collapses = collapse_if_improves_energy(_map_state, _settings);
    timer_collapse.stop();

    ISM_DEBUG_OUT("collapse done");

    ISM_ASSERT(_settings.w_approx <= 0.0 || sanity_check_all_vertices_assigned(_map_state));

    // Garbage collection
    std::vector<VH> v_map;
    std::vector<HEH> h_map;
    std::vector<FH> f_map;
    garbage_collection_with_index_maps(_map_state.mesh_T, v_map, h_map, f_map);
    for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
    {
        _map_state.maps_Tf_inputvs[i] = apply_index_map(_map_state.mesh_T, _map_state.maps_Tf_inputvs[i], f_map);
        _map_state.embeddings_T[i] = apply_index_map(_map_state.mesh_T, _map_state.embeddings_T[i], v_map);
    }

    // Update landmarks of T if landmark constraints
    for (int i = 0; i < (int)_map_state.landmarks_T.size(); i++)
    {
        auto vh_new = v_map[_map_state.landmarks_T[i].idx()];
        if (vh_new.is_valid())
            _map_state.landmarks_T[i] = vh_new;
    }

    ISM_ASSERT(_settings.w_approx <= 0.0 || sanity_check_all_vertices_assigned(_map_state));
    ISM_ASSERT(_settings.w_approx <= 0.0 || sanity_check_assignments(_map_state));

    // Release status, otherwise does not work for some reason
    _map_state.mesh_T.release_vertex_status();
    //_map_state.mesh_T.release_halfedge_status();
    _map_state.mesh_T.release_edge_status();
    _map_state.mesh_T.release_face_status();

    // Callback
    _callback("after_collapses");

    // Flip edges if improves energy
    TinyAD::Timer timer_flip("timer flip");
    if (_settings.allow_flips)
        n_flips = flip_if_improves_energy(_map_state, _settings);
    timer_flip.stop();

    ISM_ASSERT_EQ(n_flipped_volume(_map_state), 0);

    // Release status, otherwise does not work for some reason
    _map_state.mesh_T.release_vertex_status();
    //_map_state.mesh_T.release_halfedge_status();
    _map_state.mesh_T.release_edge_status();
    _map_state.mesh_T.release_face_status();

    // garbage collection
    for (int i = 0; i < (int)_map_state.meshes_input.size(); i++)
        garbage_collection(_map_state.mesh_T, _map_state.embeddings_T[i]);

    ISM_ASSERT(_settings.w_approx <= 0.0 || sanity_check_all_vertices_assigned(_map_state));
    ISM_ASSERT(_settings.w_approx <= 0.0 || sanity_check_assignments(_map_state));

    // Callback
    _callback("after_flips");

    ISM_DEBUG_OUT("Splits: " << n_splits);
    ISM_DEBUG_OUT("Collapses: "<< n_collapses);
    ISM_DEBUG_OUT("Flips: " << n_flips);
    ISM_DEBUG_OUT("#Vertices: " << _map_state.mesh_T.n_vertices());

    ISM_ASSERT(pairwise_distinct(_map_state.landmarks_T));

    // Write timings to csv
    if (!_csv_path.empty())
    {
//                "iteration", "iter_type", "derivatives (seconds)", "solve (seconds)", "line search (seconds)", "iteration (seconds)", "objective",
//                "splits (seconds)", "n_splits", "collapses (seconds)", "n_collapses", "flips (seconds)", "n_flips", "n_vertices"
        append_to_csv(_csv_path, "", "remeshing", "", "", "", "", "",
                          timer_split.seconds(), n_splits, timer_collapse.seconds(), n_collapses, timer_flip.seconds(), n_flips, _map_state.mesh_T.n_vertices());
    }

    if (n_splits + n_collapses + n_flips == 0)
        return false;

    return true;
}

}
