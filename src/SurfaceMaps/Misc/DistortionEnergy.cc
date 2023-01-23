/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt
 */

#include "DistortionEnergy.hh"

#include <TinyAD/Scalar.hh>
#include <TinyAD/Utils/Helpers.hh>
#include <SurfaceMaps/Misc/ConstantCurvatureGeometry.hh>

namespace SurfaceMaps
{

template <typename T>
T sphere_injectivity_barrier(
        const Vec3<T>& _a_sphere,
        const Vec3<T>& _b_sphere,
        const Vec3<T>& _c_sphere)
{
    // Exact injectivity check
    if (!ccw_exclusive_3d(TinyAD::to_passive(_a_sphere), TinyAD::to_passive(_b_sphere), TinyAD::to_passive(_c_sphere)))
        return INF_DOUBLE;

    // Barrier term: Oriented volume spanned by origin and triangle must be positive.
    T volume = (1.0 / 6.0) * TinyAD::col_mat(_a_sphere, _b_sphere, _c_sphere).determinant();

    // Inexact injectivity check
    if (volume <= 0.0)
        return INF_DOUBLE;

    return -log(volume);
}

// Explicit instantiation
template double sphere_injectivity_barrier(const Vec3<double>& _a_sphere, const Vec3<double>& _b_sphere, const Vec3<double>& _c_sphere);
template TinyAD::Double<6, false> sphere_injectivity_barrier( const Vec3<TinyAD::Double<6, false>>& _a_sphere, const Vec3<TinyAD::Double<6, false>>& _b_sphere, const Vec3<TinyAD::Double<6, false>>& _c_sphere);
template TinyAD::Double<6, true> sphere_injectivity_barrier(const Vec3<TinyAD::Double<6, true>>& _a_sphere, const Vec3<TinyAD::Double<6, true>>& _b_sphere, const Vec3<TinyAD::Double<6, true>>& _c_sphere);

template <typename T>
T sphere_edge_length_barrier(
        const Vec3<T>& _a_sphere,
        const Vec3<T>& _b_sphere,
        const Vec3<T>& _c_sphere,
        const double _max_edge_length_degrees)
{
    const T l_ab = distance(_a_sphere, _b_sphere, Spherical);
    const T l_bc = distance(_b_sphere, _c_sphere, Spherical);
    const T l_ca = distance(_c_sphere, _a_sphere, Spherical);

    const double threshold = _max_edge_length_degrees / 180.0 * M_PI;

    if (l_ab >= threshold || l_bc >= threshold || l_ca >= threshold)
        return INF_DOUBLE;

    ISM_ASSERT_G(threshold - l_ab, 0.0);
    ISM_ASSERT_G(threshold - l_bc, 0.0);
    ISM_ASSERT_G(threshold - l_ca, 0.0);

    return -log(threshold - l_ab) - log(threshold - l_bc) - log(threshold - l_ca);
}

// Explicit instantiation
template double sphere_edge_length_barrier(const Vec3<double>& _a_sphere, const Vec3<double>& _b_sphere, const Vec3<double>& _c_sphere, const double);
template TinyAD::Double<6, false> sphere_edge_length_barrier( const Vec3<TinyAD::Double<6, false>>& _a_sphere, const Vec3<TinyAD::Double<6, false>>& _b_sphere, const Vec3<TinyAD::Double<6, false>>& _c_sphere, const double);
template TinyAD::Double<6, true> sphere_edge_length_barrier(const Vec3<TinyAD::Double<6, true>>& _a_sphere, const Vec3<TinyAD::Double<6, true>>& _b_sphere, const Vec3<TinyAD::Double<6, true>>& _c_sphere, const double);


template <typename T>
T surface_to_sphere_angle_energy(
        const FH& _fh,
        const Vec3<T>& _a_sphere,
        const Vec3<T>& _b_sphere,
        const Vec3<T>& _c_sphere,
        const TriMesh& _mesh)
{
    HEH heh_a, heh_b, heh_c;
    handles(_mesh, _fh, heh_a, heh_b, heh_c);

    const double alpha_rest = _mesh.calc_sector_angle(heh_a);
    const double beta_rest = _mesh.calc_sector_angle(heh_b);
    const double gamma_rest = _mesh.calc_sector_angle(heh_c);

    const Vec3<T> n_ab = _b_sphere.cross(_a_sphere);
    const Vec3<T> n_bc = _c_sphere.cross(_b_sphere);
    const Vec3<T> n_ca = _a_sphere.cross(_c_sphere);

    const T alpha = atan2(n_ab.cross(-n_ca).norm(), -n_ab.dot(n_ca));
    const T beta = atan2(n_bc.cross(-n_ab).norm(), -n_bc.dot(n_ab));
    const T gamma = atan2(n_ca.cross(-n_bc).norm(), -n_ca.dot(n_bc));

    if (alpha <= 0.0 || beta <= 0.0 || gamma <= 0.0)
        return INF_DOUBLE;

    return sqr(alpha / alpha_rest) + sqr(beta / beta_rest) + sqr(gamma / gamma_rest)
         + sqr(alpha_rest / alpha) + sqr(beta_rest / beta) + sqr(gamma_rest / gamma);
}

// Explicit instantiation
template double surface_to_sphere_angle_energy(const FH& _fh, const Vec3<double>& _a_sphere, const Vec3<double>& _b_sphere, const Vec3<double>& _c_sphere, const TriMesh& _mesh);
template TinyAD::Double<6, false> surface_to_sphere_angle_energy(const FH& _fh, const Vec3<TinyAD::Double<6, false>>& _a_sphere, const Vec3<TinyAD::Double<6, false>>& _b_sphere, const Vec3<TinyAD::Double<6, false>>& _c_sphere, const TriMesh& _mesh);
template TinyAD::Double<6, true> surface_to_sphere_angle_energy(const FH& _fh, const Vec3<TinyAD::Double<6, true>>& _a_sphere, const Vec3<TinyAD::Double<6, true>>& _b_sphere, const Vec3<TinyAD::Double<6, true>>& _c_sphere, const TriMesh& _mesh);

template <typename T>
T surface_to_sphere_area_energy(
        const FH& _fh,
        const Vec3<T>& _a_sphere,
        const Vec3<T>& _b_sphere,
        const Vec3<T>& _c_sphere,
        const TriMesh& _mesh,
        const double _total_area)
{
    constexpr double sphere_area = 4.0 * M_PI;
    const double rel_area_rest = _mesh.calc_face_area(_fh) / _total_area;

    // Approximate spherical area by Euclidean area
    const T rel_area = (0.5 / sphere_area) * (_b_sphere - _a_sphere).cross(_c_sphere - _a_sphere).norm();

    ISM_ASSERT_G(rel_area, 0.0);

    return sqr(rel_area / rel_area_rest) + sqr(rel_area_rest / rel_area);
}

// Explicit instantiation
template double surface_to_sphere_area_energy(const FH& _fh, const Vec3<double>& _a_sphere, const Vec3<double>& _b_sphere, const Vec3<double>& _c_sphere, const TriMesh& _mesh, const double _total_area);
template TinyAD::Double<6, false> surface_to_sphere_area_energy(const FH& _fh, const Vec3<TinyAD::Double<6, false>>& _a_sphere, const Vec3<TinyAD::Double<6, false>>& _b_sphere, const Vec3<TinyAD::Double<6, false>>& _c_sphere, const TriMesh& _mesh, const double _total_area);
template TinyAD::Double<6, true> surface_to_sphere_area_energy(const FH& _fh, const Vec3<TinyAD::Double<6, true>>& _a_sphere, const Vec3<TinyAD::Double<6, true>>& _b_sphere, const Vec3<TinyAD::Double<6, true>>& _c_sphere, const TriMesh& _mesh, const double _total_area);

template <typename T>
T sphere_distortion_energy(
        const FH& _fh,
        const Vec3<T>& _a_sphere,
        const Vec3<T>& _b_sphere,
        const Vec3<T>& _c_sphere,
        const TriMesh& _mesh,
        const double _total_area,
        const double _w_barrier,
        const double _w_angle,
        const double _w_area,
        const double _max_edge_length_degrees)
{
    T E = 0.0;

    if (_w_barrier > 0.0)
    {
        E += _w_barrier * sphere_injectivity_barrier(_a_sphere, _b_sphere, _c_sphere)
           + _w_barrier * sphere_edge_length_barrier(_a_sphere, _b_sphere, _c_sphere, _max_edge_length_degrees);
    }

    if (_w_angle > 0.0)
        E += _w_angle * surface_to_sphere_angle_energy(_fh, _a_sphere, _b_sphere, _c_sphere, _mesh);

    if (_w_area > 0.0)
        E += _w_area * surface_to_sphere_area_energy(_fh, _a_sphere, _b_sphere, _c_sphere, _mesh, _total_area);

    return E;
}

template double sphere_distortion_energy(const FH&, const Vec3<double>&, const Vec3<double>&, const Vec3<double>&, const TriMesh&, const double, const double, const double, const double, const double);
template TinyAD::Double<2, false> sphere_distortion_energy(const FH&, const Vec3<TinyAD::Double<2, false>>&, const Vec3<TinyAD::Double<2, false>>&, const Vec3<TinyAD::Double<2, false>>&, const TriMesh&, const double, const double, const double, const double, const double);
template TinyAD::Double<2, true> sphere_distortion_energy(const FH&, const Vec3<TinyAD::Double<2, true>>&, const Vec3<TinyAD::Double<2, true>>&, const Vec3<TinyAD::Double<2, true>>&, const TriMesh&, const double, const double, const double, const double, const double);
template TinyAD::Double<6, false> sphere_distortion_energy(const FH&, const Vec3<TinyAD::Double<6, false>>&, const Vec3<TinyAD::Double<6, false>>&, const Vec3<TinyAD::Double<6, false>>&, const TriMesh&, const double, const double, const double, const double, const double);
template TinyAD::Double<6, true> sphere_distortion_energy(const FH&, const Vec3<TinyAD::Double<6, true>>&, const Vec3<TinyAD::Double<6, true>>&, const Vec3<TinyAD::Double<6, true>>&, const TriMesh&, const double, const double, const double, const double, const double);

}
