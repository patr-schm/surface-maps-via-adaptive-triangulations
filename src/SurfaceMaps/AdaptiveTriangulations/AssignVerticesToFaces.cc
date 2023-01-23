/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Dörte Pieper, Patrick Schmidt
 */

#include <SurfaceMaps/AdaptiveTriangulations/AssignVerticesToFaces.hh>

#include <SurfaceMaps/AdaptiveTriangulations/LiftToSurface.hh>
#include <SurfaceMaps/AdaptiveTriangulations/Helpers.hh>
#include <SurfaceMaps/Utils/BSPTree.hh>
#include <TinyAD/Utils/Timer.hh>
#include <queue>

namespace SurfaceMaps
{

double time_search_for_correct_face = 0.0;

SFH search_for_correct_face_for_vertex(
        const Vec3d& _p_X,      //point of mesh X on sphere
        const SFH& _start_f_X,
        const TriMesh& _mesh_embedding_X)
{
    TinyAD::Timer timer_fct("search_for_correct_face_for_vertex()", true);

    //check _start_f_T
    SVH v0, v1, v2;
    handles(_mesh_embedding_X, _start_f_X, v0, v1, v2);
    if (in_triangle_inclusive_3d(_p_X, _mesh_embedding_X.point(v0), _mesh_embedding_X.point(v1), _mesh_embedding_X.point(v2)))
    {
        //ISM_HIGHLIGHT("correct face " << _start_f_X);
        timer_fct.stop();
        time_search_for_correct_face += timer_fct.seconds();
        return _start_f_X;
    }
    //ISM_WARNING("point not in "  << _start_f_X << ", search start");
    // start face was not correct, breadth first search
    ExternalProperty<FH, bool> already_checked(_mesh_embedding_X, false);
    already_checked[_start_f_X] = true;
    std::queue<SFH> faces_to_check;
    //enqueue adjacent faces of _start_f_T
    faces_to_check.push(_start_f_X.halfedge().opp().face());
    faces_to_check.push(_start_f_X.halfedge().next().opp().face());
    faces_to_check.push(_start_f_X.halfedge().next().next().opp().face());

    while (!faces_to_check.empty())
    {
        SFH cur_face = faces_to_check.front();
        SVH v0, v1, v2;
        handles(_mesh_embedding_X, cur_face, v0, v1, v2);
        if (in_triangle_inclusive_3d(_p_X, _mesh_embedding_X.point(v0), _mesh_embedding_X.point(v1), _mesh_embedding_X.point(v2)))
        {
            timer_fct.stop();
            time_search_for_correct_face += timer_fct.seconds();
            return cur_face;
        }
        //enqueue neighbours of cur_face, if not already checked before
        auto f0 = cur_face.halfedge().opp().face();
        auto f1 = cur_face.halfedge().next().opp().face();
        auto f2 = cur_face.halfedge().next().next().opp().face();
        if (!already_checked[f0])
            faces_to_check.push(f0);
        if (!already_checked[f1])
            faces_to_check.push(f1);
        if (!already_checked[f2])
            faces_to_check.push(f2);
        already_checked[cur_face] = true;
        faces_to_check.pop();
    }
    ISM_DEBUG_OUT("NO CORRECT FACEHANDLE FOUND FOR VERTEX");
    timer_fct.stop();
    time_search_for_correct_face += timer_fct.seconds();

    return _start_f_X;
}

void assign_vertices_to_T_faces(
        MapState& _map_state)
{
    const int n_meshes = _map_state.meshes_input.size();
    _map_state.maps_Tf_inputvs.clear();
    for (int i = 0; i < n_meshes; i++)
        _map_state.maps_Tf_inputvs.push_back(ExternalProperty<FH, std::vector<SVH>> (_map_state.mesh_T));

    TinyAD::Timer timer_assign("Initial vertex assignment");
    for (int i = 0; i < n_meshes; i++)
    {
        if ((int)_map_state.active_for_approx.size() != n_meshes)
            assign_vertices_to_T_faces(_map_state.maps_Tf_inputvs[i], _map_state.meshes_input[i], _map_state.mesh_T, _map_state.embeddings_input[i], _map_state.embeddings_T[i]);
        else
            assign_vertices_to_T_faces(_map_state.maps_Tf_inputvs[i], _map_state.meshes_input[i], _map_state.mesh_T, _map_state.embeddings_input[i], _map_state.embeddings_T[i], _map_state.active_for_approx[i]);
    }

    timer_assign.stop();
    ISM_ASSERT(sanity_check_all_vertices_assigned(_map_state));
    ISM_ASSERT(sanity_check_assignments(_map_state));
}

///assign each vertex of mesh_X to the face of mesh_T based on in which it is in the sphere embedding
void assign_vertices_to_T_faces(
        ExternalProperty<FH, std::vector<SVH>>& _assign_T_X,
        const TriMesh& _mesh_X,
        const TriMesh& _mesh_T,
        const ExternalProperty<VH, Vec3d>& _embedding_X,
        const ExternalProperty<VH, Vec3d>& _embedding_T_X,
        const ExternalProperty<VH, bool>& _active_vertices_X)
{
    _assign_T_X = {_mesh_T, std::vector<SVH>()};
    TriMesh mesh_embedding_T_X = embedding_to_mesh(_mesh_T, _embedding_T_X);
    BSPTree bsp_T_X = BSPTree (mesh_embedding_T_X);

    for (auto v_X : _mesh_X.vertices())
    {
        if (_active_vertices_X.size() == _mesh_X.n_vertices())
            if (!_active_vertices_X[v_X])
                continue;
        SFH f_T = bsp_T_X.face(_embedding_X[v_X], mesh_embedding_T_X);

        // search for correct face
        f_T = search_for_correct_face_for_vertex(_embedding_X[v_X], f_T, mesh_embedding_T_X);
        _assign_T_X[f_T].push_back(v_X);
    }
}

/// Determine which of the _in_vhs lie in the triangle (on sphere), and delete them from _in_vhs
void assign_vertices_to_triangle(
        std::vector<SVH>& _in_vhs,
        std::vector<SVH>& _assigned_vhs,
        const Vec3d& _a_sphere,
        const Vec3d& _b_sphere,
        const Vec3d& _c_sphere,
        const ExternalProperty<VH, Vec3d>& _embedding_input)
{
    _assigned_vhs.clear();

    std::vector<SVH> unassigned_vhs;
    for (auto vh : _in_vhs)
    {
        if (in_triangle_inclusive_3d(_embedding_input[vh], _a_sphere, _b_sphere, _c_sphere))
            _assigned_vhs.push_back(vh);
        else
            unassigned_vhs.push_back(vh);
    }

    _in_vhs = unassigned_vhs;
}

void update_assignment_vertices_to_T_faces(
        ExternalProperty<FH, std::vector<SVH>>& _assign_T_X,
        const TriMesh& _mesh_T,
        const ExternalProperty<VH, Vec3d>& _embedding_X,
        const ExternalProperty<VH, Vec3d>& _embedding_T_X)
{
    ExternalProperty<FH, std::vector<SVH>> assign_before = _assign_T_X;
    _assign_T_X = {_mesh_T, std::vector<SVH>()};
    TriMesh mesh_embedding_T_X = embedding_to_mesh(_mesh_T, _embedding_T_X);

    for (auto fh : _mesh_T.faces())
    {
        for (auto vh : assign_before[fh])
        {
            auto new_f = search_for_correct_face_for_vertex(_embedding_X[vh], fh, mesh_embedding_T_X);
            _assign_T_X[new_f].push_back(vh);
        }
    }
}

bool sanity_check_assignments(
        const MapState& _map_state,
        const std::vector<ExternalProperty<VH, Vec3d>>& _embeddings_T)
{
    std::vector<SVH> wrong_vh_A, wrong_vh_B;
    ExternalProperty<FH, Color> _colors_T_A (_map_state.mesh_T, Color(1.0, 1.0, 1.0, 1.0));
    ExternalProperty<FH, Color> _colors_T_B (_map_state.mesh_T, Color(1.0, 1.0, 1.0, 1.0));
    for(int i = 0; i < (int)_map_state.meshes_input.size(); i++)
        ISM_ASSERT(_map_state.maps_Tf_inputvs[i].size() == _map_state.mesh_T.n_faces());

    for (auto fh : _map_state.mesh_T.faces())
    {
        SVH v0, v1, v2;
        handles(_map_state.mesh_T, fh, v0, v1, v2);
        for(int i = 0; i < (int)_map_state.meshes_input.size(); i++)
        {
            for (auto vh_i : _map_state.maps_Tf_inputvs[i][fh])
            {
                if (! in_triangle_inclusive_3d(_map_state.embeddings_input[i][vh_i], _embeddings_T[i][v0], _embeddings_T[i][v1], _embeddings_T[i][v2]))
                {
                    wrong_vh_A.push_back(vh_i);
                    ISM_DEBUG_OUT("wrong vh " << vh_i << "of mesh number " << i << " not in T face " << fh);

                    _colors_T_A[fh] = Color(0.0, 0.0, 1.0, 1.0);
                }
            }
        }
    }
    return wrong_vh_A.empty() && wrong_vh_B.empty();
}

bool sanity_check_assignments(const MapState& _map_state)
{
    return sanity_check_assignments(_map_state, _map_state.embeddings_T);
}

bool sanity_check_all_vertices_assigned(
        const MapState& _map_state)
{
    for(int i = 0; i < (int)_map_state.meshes_input.size(); i++)
    {
        int num_assigned_i = 0;
        for (auto fh : _map_state.mesh_T.faces())
        {
            num_assigned_i += _map_state.maps_Tf_inputvs[i][fh].size();
        }

        if (_map_state.active_for_approx.size() != _map_state.meshes_input.size())
        {
            if (num_assigned_i != (int)_map_state.meshes_input[i].n_vertices())
                return false;
        }
        else
        {
            return true; // For speedup just assume this is met
            // Not all vertices active so only number of active ones should be assigned
            int num_active = 0;
            for (auto vh : _map_state.meshes_input[i].vertices())
            {
                if (_map_state.active_for_approx[i][vh] == true)
                    num_active++;
            }
            if (num_assigned_i != num_active)
                return false;
        }
    }
    return true;
}

}
