/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Dörte Pieper, Patrick Schmidt
 */
#pragma once

#include <TinyAD/Utils/ToPassive.hh>
#include <SurfaceMaps/Types.hh>
#include <SurfaceMaps/Utils/Helpers.hh>
#include <SurfaceMaps/Utils/BSPTree.hh>
#include <SurfaceMaps/Misc/ConstantCurvatureGeometry.hh>
#include <SurfaceMaps/AdaptiveTriangulations/AssignVerticesToFaces.hh>

namespace SurfaceMaps
{

ExternalProperty<VH, Vec3d> lift_to_surface(
        const TriMesh& _mesh_input, // A or B
        const TriMesh& _mesh_T,
        const ExternalProperty<VH, Vec3d>& _embedding_input, // embedding_A or embedding_B
        const ExternalProperty<VH, Vec3d>& _embedding_T,
        const BSPTree& _bsp_input);

ExternalProperty<VH, Vec3d> lift_to_surface(
        const TriMesh& _mesh_input, // A or B
        const TriMesh& _mesh_T,
        const ExternalProperty<VH, Vec3d>& _embedding_input, // embedding_A or embedding_B
        const ExternalProperty<VH, Vec3d>& _embedding_T);

template <typename T>
void points(
        const TriMesh& _mesh, const FH _fh,
        Vec3<T>& _out_p_a, Vec3<T>& _out_p_b, Vec3<T>& _out_p_c)
{
    VH vh_a, vh_b, vh_c;
    handles(_mesh, _fh, vh_a, vh_b, vh_c);

    _out_p_a = _mesh.point(vh_a);
    _out_p_b = _mesh.point(vh_b);
    _out_p_c = _mesh.point(vh_c);
}

template <typename T>
void bsp_tree_barys_face(
        const Vec3<T>& _p_sphere,
        const TriMesh& _mesh_input_embedding,
        const BSPTree& _bsp_input,
        T& _alpha, T& _beta, T& _gamma,
        SFH& _f,
        const bool& _search_correct_face = true)
{
    // Find corresponding triangle of embedding on sphere, convert to Vec3d for the call
    _f = _bsp_input.face(TinyAD::to_passive(_p_sphere), _mesh_input_embedding);

    // Might be not correct face because not exact
    if (_search_correct_face)
    {
        _f = search_for_correct_face_for_vertex(TinyAD::to_passive(_p_sphere), _f, _mesh_input_embedding);
    }

    // Positions of input vertices on sphere
    Vec3<T> a_emb, b_emb, c_emb;
    points(_mesh_input_embedding, _f, a_emb, b_emb, c_emb);

    // Barycentric coordinates
    barys_abc_3d(a_emb, b_emb, c_emb, _p_sphere, Spherical, _alpha, _beta);
    _gamma =  1.0 - _alpha - _beta;
}

template <typename T>
Vec3<T> lift_vertex_to_surface(
        const Vec3<T>& _p_sphere,
        const TriMesh& _mesh_input,
        const TriMesh& _mesh_input_embedding,
        const BSPTree& _bsp_input)
{
    SFH f;
    T alpha, beta, gamma;
    bsp_tree_barys_face(_p_sphere, _mesh_input_embedding, _bsp_input, alpha, beta, gamma, f);

    // position of vertex on input mesh
    Vec3<T> a_input, b_input, c_input;
    points(_mesh_input, f, a_input, b_input, c_input);
    Vec3<T> lifted_vertex = a_input * alpha + b_input * beta + c_input * gamma;
    return lifted_vertex;
}

}
