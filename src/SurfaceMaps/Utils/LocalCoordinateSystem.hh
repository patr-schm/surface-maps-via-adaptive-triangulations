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
#include <SurfaceMaps/Utils/Helpers.hh>

namespace SurfaceMaps
{

/**
 * Local 2d coordinate system of a triangle.
 * Stores origin and basis vectors.
 * Provides convenience methods for conversion from/to the local CS.
 */
struct LocalCoordinateSystem
{

/**
 * Default constructor necessary to store this in a property.
 */
LocalCoordinateSystem()
{
    valid_ = false;
}

/**
 * Use this constructor for valid initialization.
 */
LocalCoordinateSystem(
        const FH _fh,
        const TriMesh& _mesh)
{
    origin_ = comp_origin(_mesh, _fh);
    comp_basis(_mesh, _fh, b0_, b1_);
    valid_ = true;
}

/**
 * Use this constructor for valid initialization.
 * Specify origin.
 */
LocalCoordinateSystem(
        const FH _fh,
        const Vec3d _origin,
        const TriMesh& _mesh)
{
    origin_ = _origin;
    comp_basis(_mesh, _fh, b0_, b1_);
    valid_ = true;
}

/**
 * Use this constructor for valid initialization.
 * a, b, c are triangle vertex positions (ccw).
 */
LocalCoordinateSystem(
        const Vec3d _a,
        const Vec3d _b,
        const Vec3d _c)
{
    origin_ = _a;
    comp_basis(_a, _b, _c, b0_, b1_);
    valid_ = true;
}

static LocalCoordinateSystem barycenter(
        const FH _fh,
        const TriMesh& _mesh)
{
    VH vh_a, vh_b, vh_c;
    handles(_mesh, _fh, vh_a, vh_b, vh_c);

    const Vec3d barycenter = (_mesh.point(vh_a) + _mesh.point(vh_b) + _mesh.point(vh_c)) / 3.0;
    return LocalCoordinateSystem(_fh, barycenter, _mesh);
}

/**
 * Returns origin of the coordinate system in model space (GCS)
 */
Vec3d origin() const
{
    return origin_;
}

/**
 * Returns the origin of the coordinate system, i.e. the vertex the first
 * half-edge is pointing from.
 */
static Vec3d comp_origin(
        const TriMesh& _mesh,
        const FH _fh)
{
    const HEH heh0 = _mesh.halfedge_handle(_fh);
    const VH vh = _mesh.from_vertex_handle(heh0);

    return _mesh.point(vh);
}

/**
 * Computes the orthonormal basis vectors of a 2d coordinate system in
 * which b0 corresponds to the direction of the first half-edge and
 * b1 is obtained by a 90° ccw rotation.
 */
static Mat<3, 2, double> comp_linear_basis(
        const TriMesh& _mesh,
        const FH& _fh)
{
    const auto normal = _mesh.calc_face_normal(_fh);
    const auto heh0 = _mesh.halfedge_handle(_fh);

    Mat<3, 2, double> B;
    B.col(0) = _mesh.calc_edge_vector(heh0).normalized();
    B.col(1) = normal.cross(B.col(0));

    ISM_ASSERT(fabs(normal.squaredNorm() - 1.0) < 1e-6);
    ISM_ASSERT(fabs(B.col(1).squaredNorm() - 1.0) < 1e-6);

    return B;
}

/**
 * Computes the orthonormal basis vectors of a 2d coordinate system in
 * which b0 corresponds to the direction of the first half-edge and
 * b1 is obtained by a 90° ccw rotation.
 */
static void comp_basis(
        const TriMesh &_mesh,
        const FH _fh,
        Vec3d& _out_b0,
        Vec3d& _out_b1)
{
    const auto normal = _mesh.calc_face_normal(_fh);
    const auto heh0 = _mesh.halfedge_handle(_fh);

    _out_b0 = _mesh.calc_edge_vector(heh0).normalized();
    _out_b1 = normal.cross(_out_b0);

    ISM_ASSERT(fabs(normal.squaredNorm() - 1.0) < 1e-6);
    ISM_ASSERT(fabs(_out_b1.squaredNorm() - 1.0) < 1e-6);
}

/**
 * Computes the orthonormal basis vectors of a 2d coordinate system in
 * which b0 corresponds to the direction of b-a and b1 is obtained
 * by a 90° ccw rotation.
 */
static void comp_basis(
        const Vec3d _a,
        const Vec3d _b,
        const Vec3d _c,
        Vec3d &_out_b0,
        Vec3d &_out_b1)
{
    const auto normal = ((_b - _a).cross(_c - _a)).normalized();

    _out_b0 = (_b - _a).normalized();
    _out_b1 = normal.cross(_out_b0);
}

/**
 * Projects a vector from model space to the local coordinate system.
 */
Vec2d vector_to_local(
        const Vec3d& _vector) const
{
    return Vec2d(_vector.dot(b0_), _vector.dot(b1_));
}

/**
 * Maps a vector from the local coordinate system to model space.
 */
Vec3d vector_to_global(
        const Vec2d& _vector_local) const
{
    return _vector_local[0] * b0_ + _vector_local[1] * b1_;
}

/**
 * Maps a point from model space to the local coordinate system.
 */
Vec2d point_to_local(
        Vec3d _point) const
{
    _point -= origin_;

    return Vec2d(_point.dot(b0_), _point.dot(b1_));
}

/**
 * Maps 3 points from model space to the local coordinate system.
 */
void pointsToLocal(
        const Vec3d _pa,
        const Vec3d _pb,
        const Vec3d _pc,
        Vec2d& _out_a,
        Vec2d& _out_b,
        Vec2d& _out_c) const
{
    _out_a = point_to_local(_pa);
    _out_b = point_to_local(_pb);
    _out_c = point_to_local(_pc);
}

/**
 * Maps a point from the local coordinate system to model space.
 */
Vec3d point_to_global(
        const Vec2d &_point_local) const
{
    return origin_ + _point_local[0] * b0_ + _point_local[1] * b1_;
}

template <typename T>
Vec3<T> point_to_global(
        const Vec2<T>& _point_local) const
{
    return origin_ + _point_local[0] * b0_ + _point_local[1] * b1_;
}

/**
 * Maps the vertices of a triangle from model space to its local coordinate system
 */
static void vertices_to_local(
        FH _fh,
        const TriMesh& _mesh,
        Vec2d& _out_a,
        Vec2d& _out_b,
        Vec2d& _out_c)
{
    VH vh_a, vh_b, vh_c;
    handles(_mesh, _fh, vh_a, vh_b, vh_c);

    const auto a_3d = _mesh.point(vh_a);
    const auto b_3d = _mesh.point(vh_b);
    const auto c_3d = _mesh.point(vh_c);

    const LocalCoordinateSystem cs(_fh, _mesh);
    _out_a = cs.point_to_local(a_3d);
    _out_b = cs.point_to_local(b_3d);
    _out_c = cs.point_to_local(c_3d);
}

Vec3d origin_;
Vec3d b0_;
Vec3d b1_;

bool valid_;

};

}
