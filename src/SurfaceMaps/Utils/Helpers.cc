/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt, Janis Born
 */

#include "Helpers.hh"

namespace SurfaceMaps
{

Vec3d center_of_gravity(
        const TriMesh& _mesh)
{
    Vec3d cog = Vec3d::Zero();
    for (auto v : _mesh.vertices())
        cog += _mesh.point(v);
    cog /= _mesh.n_vertices();

    return cog;
}

Vec3d center_of_gravity(
        const TriMesh& _mesh,
        const ExternalProperty<VH, Vec3d>& _embedding)
{
    Vec3d cog = Vec3d::Zero();
    for (auto v : _mesh.vertices())
        cog += _embedding[v];
    cog /= _mesh.n_vertices();

    return cog;
}

double total_area(
        const TriMesh &_mesh)
{
    double result = 0.0;
    for (auto fh : _mesh.faces())
        result += _mesh.calc_face_area(fh);

    return result;
}

double total_area(
        const TriMesh& _mesh,
        const ExternalProperty<VH, Vec3d>& _embedding)
{
    double result = 0.0;
    for (auto fh : _mesh.faces())
    {
        VH vh_a, vh_b, vh_c;
        handles(_mesh, fh, vh_a, vh_b, vh_c);
        const Vec3d a = _embedding[vh_a];
        const Vec3d b = _embedding[vh_b];
        const Vec3d c = _embedding[vh_c];

        result += 0.5 * (b-a).cross(c-a).norm();
    }

    return result;
}

ExternalProperty<VH, double> vertex_areas(
        const TriMesh& _mesh)
{
    ExternalProperty<VH, double> areas(_mesh, 0.0);
    for (auto f : _mesh.faces())
    {
        const double area = _mesh.calc_face_area(f);
        for (auto v : f.vertices())
            areas[v] += area / 3.0;
    }

    return areas;
}

ExternalProperty<VH, double> rel_vertex_areas(
        const TriMesh& _mesh)
{
    const double total = total_area(_mesh);

    ExternalProperty<VH, double> areas(_mesh, 0.0);
    for (auto f : _mesh.faces())
    {
        const double area = _mesh.calc_face_area(f);
        for (auto v : f.vertices())
            areas[v] += area / 3.0 / total;
    }

    return areas;
}

std::pair<Vec3d, Vec3d> bounding_box(
        const TriMesh& _mesh)
{
    Vec3d pmin(+INF_DOUBLE, +INF_DOUBLE, +INF_DOUBLE);
    Vec3d pmax(-INF_DOUBLE, -INF_DOUBLE, -INF_DOUBLE);
    for (const auto& vh : _mesh.vertices())
    {
        const Vec3d& p = _mesh.point(vh);
        pmin = pmin.cwiseMin(p);
        pmax = pmax.cwiseMax(p);
    }
    return {pmin, pmax};
}

double bounding_box_diagonal(
        const TriMesh& _mesh)
{
    const auto [pmin, pmax] = bounding_box(_mesh);
    return (pmax - pmin).norm();
}

bool incident(const TriMesh& _mesh, const VH _vh, const FH _fh)
{
    ISM_ASSERT(_mesh.is_valid_handle(_fh));
    for (const auto& vh : _mesh.fv_range(_fh))
        if (vh == _vh)
            return true;
    return false;
}

VH closest_vertex(
        const TriMesh& _mesh,
        const Vec3d& _p)
{
    double closest_dist = INF_DOUBLE;
    VH closest_vh;
    for (auto v : _mesh.vertices())
    {
        const double dist = (_mesh.point(v) - _p).norm();
        if (dist < closest_dist)
        {
            closest_dist = dist;
            closest_vh = v;
        }
    }

    return closest_vh;
}

void rotate(
        TriMesh& _mesh,
        const double _angle,
        const Vec3d& _axis)
{
     Mat3d R = Eigen::AngleAxis(_angle, _axis.normalized()).matrix();
     for (auto v : _mesh.vertices())
         _mesh.point(v) = R * _mesh.point(v);
}

void rotate(
        const TriMesh& _mesh,
        ExternalProperty<VH, Vec3d>& _embedding,
        const double _angle,
        const Vec3d& _axis)
{
    Mat3d R = Eigen::AngleAxis(_angle, _axis.normalized()).matrix();
    for (auto v : _mesh.vertices())
        _embedding[v] = R * _embedding[v];
}

void translate(
        TriMesh& _mesh,
        const Vec3d& _t)
{
     for (auto v : _mesh.vertices())
         _mesh.point(v) += _t;
}

}
