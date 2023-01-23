/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt
 */
#pragma once

#include <SurfaceMaps/Types.hh>
#include <SurfaceMaps/Utils/Geometries.hh>
#include <SurfaceMaps/Utils/Helpers.hh>

namespace SurfaceMaps
{

template <typename T>
T distance(const Vec3<T>& _a, const Vec3<T>& _b, const Geometry& _geometry)
{
    if (_geometry == Spherical)
        return atan2(_a.cross(_b).norm(), _a.dot(_b));
    else
        ISM_ERROR_throw("");
}

template <typename T>
Vec3<T> project_to_model(const Vec3<T>& _p, const Geometry& _geometry)
{
    if (_geometry == Spherical)
        return _p.normalized();
    else
        ISM_ERROR_throw("");
}

// Compute barycentric coordinates in projected triangle a', b', c'
// at intersection point with p.
template <typename T1, typename T2, typename TR>
void barys_abc_3d(
        const Vec3<T1>& a, const Vec3<T1>& b, const Vec3<T1>& c, // Do not have to lie on manifold
        const Vec3<T2>& p, // Does lie on manifold
        const Geometry& _geometry,
        TR& alpha,
        TR& beta)
{
    ISM_ASSERT(_geometry == Planar || _geometry == Spherical || _geometry == Hyperbolic);

    const Vec3<T1> a_pr = project_to_model(a, _geometry);
    const Vec3<T1> b_pr = project_to_model(b, _geometry);
    const Vec3<T1> c_pr = project_to_model(c, _geometry);

    Mat3<TR> M_alpha;
    M_alpha << -c_pr, b_pr - c_pr, -p;

    Mat3<TR> M_beta;
    M_beta << a_pr - c_pr, - c_pr, -p;

    Mat3<TR> M;
    M << a_pr - c_pr, b_pr - c_pr, -p;

    TR M_det = M.determinant();

    // Cramer's rule for 3x3 system
    alpha = M_alpha.determinant() / M_det;
    beta = M_beta.determinant() / M_det;
}

bool ccw_inclusive_3d(
        const Vec3d& a, const Vec3d& b, const Vec3d& c);

bool ccw_exclusive_3d(
        const Vec3d& a, const Vec3d& b, const Vec3d& c);

bool in_triangle_inclusive_3d(
        const Vec3d& p,
        const Vec3d& a, const Vec3d& b, const Vec3d& c);

bool in_triangle_exclusive_3d(
        const Vec3d& p,
        const Vec3d& a, const Vec3d& b, const Vec3d& c);

}
