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

namespace SurfaceMaps
{

template <typename T>
T sphere_injectivity_barrier(
        const Vec3<T>& _a_sphere,
        const Vec3<T>& _b_sphere,
        const Vec3<T>& _c_sphere);

template <typename T>
T sphere_edge_length_barrier(
        const Vec3<T>& _a_sphere,
        const Vec3<T>& _b_sphere,
        const Vec3<T>& _c_sphere,
        const double _max_edge_length_degrees);

template <typename T>
T surface_to_sphere_angle_energy(
        const FH& _fh,
        const Vec3<T>& _a_sphere,
        const Vec3<T>& _b_sphere,
        const Vec3<T>& _c_sphere,
        const TriMesh& _mesh);

template <typename T>
T surface_to_sphere_area_energy(
        const FH& _fh,
        const Vec3<T>& _a_sphere,
        const Vec3<T>& _b_sphere,
        const Vec3<T>& _c_sphere,
        const TriMesh& _mesh,
        const double _total_area);

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
        const double _max_edge_length_degrees);

}
