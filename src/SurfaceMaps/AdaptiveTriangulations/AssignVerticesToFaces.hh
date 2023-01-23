/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Dörte Pieper, Patrick Schmidt
 */
#pragma once

#include <SurfaceMaps/Types.hh>
#include <SurfaceMaps/AdaptiveTriangulations/MapState.hh>
#include <SurfaceMaps/AdaptiveTriangulations/AdaptiveTriangulationsSettings.hh>

namespace SurfaceMaps
{

extern double time_search_for_correct_face;

/// breadth first search for face in which vertex is
SFH search_for_correct_face_for_vertex(
        const Vec3d& _p_X,      //point of mesh X on sphere
        const SFH& _start_f_T,
        const TriMesh& _mesh_embedding_T_X);

void assign_vertices_to_T_faces(
        MapState& _map_state);

///assign each vertice of mesh_X to the face of mesh_T based on in which it is in the sphere embedding
void assign_vertices_to_T_faces(
        ExternalProperty<FH, std::vector<SVH>>& _assign_T_X,
        const TriMesh& _mesh_X,
        const TriMesh& _mesh_T,
        const ExternalProperty<VH, Vec3d>& _embedding_X,
        const ExternalProperty<VH, Vec3d>& _embedding_T_X,
        const ExternalProperty<VH, bool>& _active_vertices_X = ExternalProperty<VH, bool>());

/// Determine which of the _in_vhs lie in the triangle (on sphere), and delete them from _in_vhs
void assign_vertices_to_triangle(
        std::vector<SVH>& _in_vhs,
        std::vector<SVH>& _assigned_vhs,
        const Vec3d& _a_sphere,
        const Vec3d& _b_sphere,
        const Vec3d& _c_sphere,
        const ExternalProperty<VH, Vec3d>& _embedding_input);

/// update assignment with local search starting from a previous assignment
void update_assignment_vertices_to_T_faces(
        ExternalProperty<FH, std::vector<SVH>>& _assign_T_X,
        const TriMesh& _mesh_T,
        const ExternalProperty<VH, Vec3d>& _embedding_X,
        const ExternalProperty<VH, Vec3d>& _embedding_T_X);

/// check MapState if for each face in mesh_T all assigned vertices in maps_Tf_inputvs are really in the face of given _embeddings_T
bool sanity_check_assignments(
        const MapState& _map_state,
        const std::vector<ExternalProperty<VH, Vec3d>>& _embeddings_T);

/// check if for each face in mesh_T all assigned vertices in map_Tf_Xvs of mesh A and mesh B are really in the face
bool sanity_check_assignments(
        const MapState& _map_state);

/// check if number of assigned vertices equals to mesh vertices (does not check if each only once)
bool sanity_check_all_vertices_assigned(
        const MapState& _map_state);

}
