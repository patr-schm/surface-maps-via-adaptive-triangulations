/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt, Janis Born
 */
#pragma once

#include <SurfaceMaps/Types.hh>

namespace SurfaceMaps
{

void center_mesh(
        TriMesh& _mesh,
        const TriMesh::Point& _new_origin = {0,0,0});

/**
 * Returns factor by which mesh was scaled
 */
double normalize_surface_area(
        TriMesh& _mesh,
        double _new_area = 1.0);

/**
 * Returns factor by which mesh was scaled
 */
double normalize_mesh(
        TriMesh& _mesh);

/**
 * Compute rigid transformation (+scaling)
 * that aligns points V_B to points V_A
 */
Eigen::Affine3d compute_rigid_alignment(
        const MatXd &V_A, // points of A as rows
        const MatXd &V_B, // points of B as rows
        bool _allow_scaling = false);

/**
 * Compute rigid transformation (+scaling)
 * that aligns mesh_B to mesh_A.
 */
Eigen::Affine3d compute_rigid_alignment(
        const std::vector<VH> &_vertices_A,
        const std::vector<VH> &_vertices_B,
        const TriMesh &_mesh_A,
        const TriMesh &_mesh_B,
        bool _allow_scaling = false);

void align_rigid(
        const std::vector<VH>& _vertices_A,
        const std::vector<VH>& _vertices_B,
        const TriMesh& _mesh_A, // reference
        TriMesh& _mesh_B, // output
        bool _allow_scaling = false);

Eigen::Affine3d compute_rotation_alignment(
        const MatXd &V_A, // points of A as rows
        const MatXd &V_B); // points of B as rows

void align_rotation(
        const std::vector<VH> &_vertices_A,
        const std::vector<VH> &_vertices_B,
        const TriMesh &_mesh_A,
        const TriMesh &_mesh_B,
        const ExternalProperty<VH, Vec3d>& _embedding_A,
        ExternalProperty<VH, Vec3d>& _embedding_B);

void normalize_and_align(
        TriMesh& _mesh_A,
        TriMesh& _mesh_B,
        const std::vector<VH>& _landmarks_A,
        const std::vector<VH>& _landmarks_B,
        const bool _normalize,
        const bool _align);

/**
 * Reflect vertex positions along one coordinate axis
 * and flip face winding order.
 */
void flip_mesh(
        TriMesh& _mesh);

}
