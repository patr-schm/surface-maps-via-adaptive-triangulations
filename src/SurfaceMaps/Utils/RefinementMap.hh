/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Janis Born
 */
#pragma once

#include <SurfaceMaps/Types.hh>
#include <SurfaceMaps/Utils/BarycentricPoint.hh>

namespace SurfaceMaps
{

/// Given some mesh M2 which is a refinement / subdivision of some mesh M1, a
/// RefinementMap assigns each face corner of M2 the corresponding original
/// location on M1.
/// Can be used to transfer data (e.g. uv coordinates) from M1 to M2 via
/// barycentric interpolation.
/// Condition: For every face of M2, all three corners of the RefinementMap
/// must map to BarycentricPoints on the same face of M1.
using RefinementMap = ExternalProperty<HEH, BarycentricPoint>;

inline
RefinementMap
identity_refinement_map(const TriMesh& _mesh)
{
    RefinementMap result(_mesh);
    for (auto f : _mesh.faces())
    {
        SHEH h_ab;
        SHEH h_bc;
        SHEH h_ca;
        handles(_mesh, f, h_ca, h_ab, h_bc);

        const BarycentricPoint bary_a(f, 1.0, 0.0, _mesh);
        const BarycentricPoint bary_b(f, 0.0, 1.0, _mesh);
        const BarycentricPoint bary_c(f, 0.0, 0.0, _mesh);

        result[h_ca] = bary_a;
        result[h_ab] = bary_b;
        result[h_bc] = bary_c;
    }
    return result;
}

template <typename T>
ExternalProperty<HEH, T>
refine_property(
        const TriMesh& _source_mesh,
        const TriMesh& _target_mesh,
        const RefinementMap& _ref_map,
        const ExternalProperty<HEH, T>& _source_prop)
{
    ExternalProperty<HEH, T> result(_target_mesh);
    for (const auto heh : _target_mesh.halfedges())
    {
        const auto& interpolated = _ref_map[heh].interpolate(_source_prop, _source_mesh);
        result[heh] = interpolated;
    }
    return result;
}

template <typename T>
ExternalProperty<VH, T>
refine_property(
        const TriMesh& _source_mesh,
        const TriMesh& _target_mesh,
        const RefinementMap& _ref_map,
        const ExternalProperty<VH, T>& _source_prop)
{
    ExternalProperty<VH, T> result(_target_mesh);
    for (const auto vh : _target_mesh.vertices())
    {
        const auto heh = vh.in();
        const auto& interpolated = _ref_map[heh].interpolate(_source_prop, _source_mesh);
        result[vh] = interpolated;
    }
    return result;
}

template <typename T>
void
refine_property(
        const TriMesh& _source_mesh,
        TriMesh& _target_mesh,
        const RefinementMap& _ref_map,
        const OpenMesh::HPropHandleT<T>& _source_prop,
        const OpenMesh::HPropHandleT<T>& _target_prop)
{
    ExternalProperty<HEH, T> ext_orig_prop(_source_mesh, _source_prop);
    ExternalProperty<HEH, T> ext_new_prop = refine_property(_source_mesh, _target_mesh, _ref_map, ext_orig_prop);
    for (const auto& heh : _target_mesh.halfedges())
    {
        _target_mesh.property(_target_prop, heh) = ext_new_prop[heh];
    }
}

template <typename T>
void
refine_property(
        const TriMesh& _source_mesh,
        TriMesh& _target_mesh,
        const RefinementMap& _ref_map,
        const OpenMesh::VPropHandleT<T>& _source_prop,
        const OpenMesh::VPropHandleT<T>& _target_prop)
{
    ExternalProperty<VH, T> ext_orig_prop(_source_mesh, _source_prop);
    ExternalProperty<VH, T> ext_new_prop = refine_property(_source_mesh, _target_mesh, _ref_map, ext_orig_prop);
    for (const auto& vh : _target_mesh.vertices())
    {
        _target_mesh.property(_target_prop, vh) = ext_new_prop[vh];
    }
}

inline
void
refine_tex_coords(
        const TriMesh& _source_mesh,
        TriMesh& _target_mesh,
        const RefinementMap& _ref_map)
{
    if (_source_mesh.has_halfedge_texcoords2D())
    {
        _target_mesh.request_halfedge_texcoords2D();
        refine_property(_source_mesh, _target_mesh, _ref_map, _source_mesh.halfedge_texcoords2D_pph(), _target_mesh.halfedge_texcoords2D_pph());
    }
    if (_source_mesh.has_vertex_texcoords2D())
    {
        _target_mesh.request_vertex_texcoords2D();
        refine_property(_source_mesh, _target_mesh, _ref_map, _source_mesh.vertex_texcoords2D_pph(), _target_mesh.vertex_texcoords2D_pph());
    }
}

}
