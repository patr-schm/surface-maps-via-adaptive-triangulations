/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt, Janis Born
 */

#include "MeshNormalization.hh"

#include <SurfaceMaps/Utils/Helpers.hh>
#include <SurfaceMaps/Utils/Genus.hh>

namespace SurfaceMaps
{

void center_mesh(
        TriMesh &_mesh,
        const TriMesh::Point &_new_origin)
{
    TriMesh::Point cog = {0,0,0};
    for (const VH vh : _mesh.vertices())
        cog += _mesh.point(vh);
    cog /= _mesh.n_vertices();
    const TriMesh::Point shift = _new_origin - cog;
    for (const VH vh : _mesh.vertices())
    {
        _mesh.point(vh) += shift;
    }

//    ISM_DEBUG_OUT(__FUNCTION__ << " moved mesh by " << shift);
}

double normalize_surface_area(
        TriMesh &_mesh,
        double _new_area)
{
    const double area = total_area(_mesh);
    const double scale = std::sqrt(_new_area) / std::sqrt(area);
    for (const VH vh : _mesh.vertices())
        _mesh.point(vh) *= scale;
    ISM_ASSERT_EPS(total_area(_mesh), _new_area, 1e-6);

//    ISM_DEBUG_OUT(__FUNCTION__ << " scaled mesh by " << scale);
    return scale;
}

double normalize_mesh(
        TriMesh& _mesh)
{
    center_mesh(_mesh, {0, 0, 0});
    return normalize_surface_area(_mesh, 1.0);
}

Eigen::Affine3d compute_rigid_alignment(
        const MatXd &V_A, // points of A as rows
        const MatXd &V_B, // points of B as rows
        bool _allow_scaling)
{
    ISM_ASSERT_EQ(V_A.cols(), 3);
    ISM_ASSERT_EQ(V_B.cols(), 3);
    ISM_ASSERT_EQ(V_A.rows(), V_B.rows());

    Eigen::Vector3d cog_A = V_A.colwise().mean();
    Eigen::MatrixXd V_A_centered = V_A;
    V_A_centered.rowwise() -= cog_A.transpose();

    Eigen::Vector3d cog_B = V_B.colwise().mean();
    Eigen::MatrixXd V_B_centered = V_B;
    V_B_centered.rowwise() -= cog_B.transpose();

    Eigen::Matrix3d cov = V_A_centered.transpose() * V_B_centered;
    Eigen::JacobiSVD<decltype(cov)> svd(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d rot = svd.matrixU() * svd.matrixV().transpose();
    if (rot.determinant() < 0.0)
    {
        Eigen::Vector3d diag(1.0, 1.0, -1.0);
        rot = svd.matrixU() * diag.asDiagonal() * svd.matrixV().transpose();
    }

    double scale = 1.0;
    if (_allow_scaling)
        scale = std::sqrt(V_A_centered.squaredNorm()) / std::sqrt(V_B_centered.squaredNorm());

    return Eigen::Translation3d(cog_A) * (scale * rot * Eigen::Translation3d(-cog_B));
}

Eigen::Affine3d compute_rigid_alignment(
        const std::vector<VH> &_vertices_A,
        const std::vector<VH> &_vertices_B,
        const TriMesh &_mesh_A,
        const TriMesh &_mesh_B,
        bool _allow_scaling)
{
    ISM_ASSERT_GEQ(_vertices_A.size(), 2);
    ISM_ASSERT_GEQ(_vertices_B.size(), 2);
    ISM_ASSERT_EQ(_vertices_A.size(), _vertices_B.size());

    Eigen::MatrixXd V_A(_vertices_A.size(), 3);
    Eigen::MatrixXd V_B(_vertices_B.size(), 3);
    for (int i = 0; i < (int)_vertices_A.size(); ++i)
    {
        V_A.row(i) = _mesh_A.point(_vertices_A[i]);
        V_B.row(i) = _mesh_B.point(_vertices_B[i]);
    }

    return compute_rigid_alignment(V_A, V_B, _allow_scaling);
}

void align_rigid(
        const std::vector<VH> &_vertices_A,
        const std::vector<VH> &_vertices_B,
        const TriMesh &_mesh_A,
        TriMesh &_mesh_B,
        bool _allow_scaling)
{
    Eigen::Affine3d trans = compute_rigid_alignment(_vertices_A, _vertices_B, _mesh_A, _mesh_B, _allow_scaling);

    for (const auto& vh : _mesh_B.vertices())
        _mesh_B.point(vh) = trans * _mesh_B.point(vh);

    _mesh_B.update_normals();
}

Eigen::Affine3d compute_rotation_alignment(
        const MatXd &V_A, // points of A as rows
        const MatXd &V_B) // points of B as rows
{
    ISM_ASSERT_EQ(V_A.cols(), 3);
    ISM_ASSERT_EQ(V_B.cols(), 3);
    ISM_ASSERT_EQ(V_A.rows(), V_B.rows());

    Eigen::Matrix3d cov = V_A.transpose() * V_B;
    Eigen::JacobiSVD<decltype(cov)> svd(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d rot = svd.matrixU() * svd.matrixV().transpose();
    if (rot.determinant() < 0.0)
    {
        Eigen::Vector3d diag(1.0, 1.0, -1.0);
        rot = svd.matrixU() * diag.asDiagonal() * svd.matrixV().transpose();
    }

    return Eigen::Affine3d(rot);
}

Eigen::Affine3d compute_rotation_alignment(
        const std::vector<VH> &_vertices_A,
        const std::vector<VH> &_vertices_B,
        const TriMesh &_mesh_A,
        const TriMesh &_mesh_B)
{
    ISM_ASSERT_GEQ(_vertices_A.size(), 4);
    ISM_ASSERT_GEQ(_vertices_B.size(), 4);
    ISM_ASSERT_EQ(_vertices_A.size(), _vertices_B.size());

    Eigen::MatrixXd V_A(_vertices_A.size(), 3);
    Eigen::MatrixXd V_B(_vertices_B.size(), 3);
    for (int i = 0; i < (int)_vertices_A.size(); ++i)
    {
        V_A.row(i) = _mesh_A.point(_vertices_A[i]);
        V_B.row(i) = _mesh_B.point(_vertices_B[i]);
    }

    return compute_rotation_alignment(V_A, V_B);
}

void align_rotation(
        const std::vector<VH> &_vertices_A,
        const std::vector<VH> &_vertices_B,
        const TriMesh &_mesh_A,
        const TriMesh &_mesh_B,
        const ExternalProperty<VH, Vec3d>& _embedding_A,
        ExternalProperty<VH, Vec3d>& _embedding_B)
{
    ISM_ASSERT_GEQ(_vertices_A.size(), 2);
    ISM_ASSERT_GEQ(_vertices_B.size(), 2);
    ISM_ASSERT_EQ(_vertices_A.size(), _vertices_B.size());

    Eigen::MatrixXd V_A(_vertices_A.size(), 3);
    Eigen::MatrixXd V_B(_vertices_B.size(), 3);
    for (int i = 0; i < (int)_vertices_A.size(); ++i)
    {
        V_A.row(i) = _embedding_A[_vertices_A[i]];
        V_B.row(i) = _embedding_B[_vertices_B[i]];
    }

    Eigen::Affine3d trans = compute_rotation_alignment(V_A, V_B);

    for (const auto& vh : _mesh_B.vertices())
        _embedding_B[vh] = trans * _embedding_B[vh];
}

namespace
{

void warn_thin_triangles(
        TriMesh& _mesh,
        const std::string& _mesh_name)
{
    // Warn if the input mesh contains thin triangles
    for (auto f : _mesh.faces())
    {
        // Find longest edge
        auto h_max = f.halfedge();
        for (auto h : f.halfedges())
        {
            if (_mesh.calc_edge_length(h) > _mesh.calc_edge_length(h_max))
                h_max = h;
        }

        // Compute aspect ratio
        const Vec3d a = _mesh.point(h_max.from());
        const Vec3d b = _mesh.point(h_max.to());
        const Vec3d c = _mesh.point(h_max.next().to());

        const Vec3d n = (b - a).normalized();
        const Vec3d p = a + n.dot(c - a) * n;

        const double base = (b - a).norm();
        const double height = (c - p).norm();

        const double ratio = base / height;

        const double aspect_ratio_thresh = 1000;
        if (ratio > aspect_ratio_thresh)
            ISM_WARNING("Input mesh " << _mesh_name << " contains thin triangle with aspect ratio " << ratio);
    }
}

}

void normalize_and_align(
        TriMesh& _mesh_A,
        TriMesh& _mesh_B,
        const std::vector<VH>& _landmarks_A,
        const std::vector<VH>& _landmarks_B,
        const bool _normalize,
        const bool _align)
{
    ISM_ASSERT(closed_surface(_mesh_A));
    ISM_ASSERT(closed_surface(_mesh_B));
    ISM_INFO("Genus A: " << genus(_mesh_A));
    ISM_INFO("Genus B: " << genus(_mesh_B));
    ISM_ASSERT_EQ(genus(_mesh_A), genus(_mesh_B));
    ISM_ASSERT_EQ(_landmarks_A.size(), _landmarks_B.size());
    ISM_ASSERT(pairwise_distinct(_landmarks_A));
    ISM_ASSERT(pairwise_distinct(_landmarks_B));

    if (_normalize)
    {
        center_mesh(_mesh_A);
        center_mesh(_mesh_B);
        normalize_surface_area(_mesh_A);
        normalize_surface_area(_mesh_B);
    }

    if (_align)
    {
        if (_landmarks_A.size() < 4)
            ISM_DEBUG_OUT("Cannot align meshes: Too few landmarks.")
        else
            align_rigid(_landmarks_B, _landmarks_A, _mesh_B, _mesh_A);
    }

    // Warn if input meshes have thin triangles
    warn_thin_triangles(_mesh_A, "A");
    warn_thin_triangles(_mesh_B, "B");

    // Check for degeneracies
    for (auto e : _mesh_A.edges())
        ISM_ASSERT_G(_mesh_A.calc_edge_length(e), 0.0);
    for (auto e : _mesh_B.edges())
        ISM_ASSERT_G(_mesh_B.calc_edge_length(e), 0.0);
}

void flip_mesh(
        TriMesh& _mesh)
{
    TriMesh mesh_flipped;

    auto reflect = [] (const Vec3d& p)
    {
        return Vec3d(p[0], p[1], -p[2]);
    };

    // Copy vertices
    for (auto v : _mesh.vertices())
        mesh_flipped.add_vertex(reflect(_mesh.point(v)));

    // Copy faces with inverted winding order
    for (auto f : _mesh.faces())
    {
        VH va, vb, vc;
        handles(_mesh, f, va, vb, vc);
        mesh_flipped.add_face(vc, vb, va);
    }

    _mesh = mesh_flipped;
}

}
