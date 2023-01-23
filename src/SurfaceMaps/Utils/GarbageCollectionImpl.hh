/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt
 */
#pragma once

#include <SurfaceMaps/Utils/GarbageCollection.hh>
#include <SurfaceMaps/Utils/ExternalProperty.hh>

namespace SurfaceMaps
{

template <typename Mesh>
void garbage_collection_with_index_maps(
        Mesh& _mesh,
        std::vector<VH>& _vertex_map,
        std::vector<HEH>& _halfedge_map,
        std::vector<FH>& _face_map)
{
    // Init identity maps
    _vertex_map.resize(_mesh.n_vertices());
    for (int i = 0; i < (int)_vertex_map.size(); ++i)
        _vertex_map[i] = VH(i);
    _halfedge_map.resize(_mesh.n_halfedges());
    for (int i = 0; i < (int)_halfedge_map.size(); ++i)
        _halfedge_map[i] = HEH(i);
    _face_map.resize(_mesh.n_faces());
    for (int i = 0; i < (int)_face_map.size(); ++i)
        _face_map[i] = FH(i);

    // Init vectors of pointers to handles
    std::vector<VH*> vh_ptrs(_vertex_map.size());
    std::vector<HEH*> heh_ptrs(_halfedge_map.size());
    std::vector<FH*> fh_ptrs(_face_map.size());
    for (int i = 0; i < (int)_vertex_map.size(); ++i)
        vh_ptrs[i] = &_vertex_map[i];
    for (int i = 0; i < (int)_halfedge_map.size(); ++i)
        heh_ptrs[i] = &_halfedge_map[i];
    for (int i = 0; i < (int)_face_map.size(); ++i)
        fh_ptrs[i] = &_face_map[i];

    // Perform garbage collection and adjust vertex maps
    // Deleted handles will be set to -1
    _mesh.garbage_collection(vh_ptrs, heh_ptrs, fh_ptrs);
}

template <typename Mesh, typename HandleType, typename ValueType>
ExternalProperty<HandleType, ValueType> apply_index_map(
        const Mesh& _mesh,
        const ExternalProperty<HandleType, ValueType>& _prop_old,
        const std::vector<HandleType>& _idx_map)
{
    ExternalProperty<HandleType, ValueType> prop_new(_mesh);
    for (int i = 0; i < (int)_prop_old.size(); ++i)
    {
        HandleType handle_new = _idx_map[i];
        if (handle_new.is_valid())
            prop_new[handle_new] = _prop_old[HandleType(i)];
    }

    return prop_new;
}

template <typename Mesh, typename ValueType>
void garbage_collection(
        Mesh& _mesh,
        ExternalProperty<VH, ValueType>& _v_prop)
{
    ISM_ASSERT(_v_prop.size_okay(_mesh));

    std::vector<VH> v_map;
    std::vector<HEH> h_map;
    std::vector<FH> f_map;
    garbage_collection_with_index_maps(_mesh, v_map, h_map, f_map);

    _v_prop = apply_index_map(_mesh, _v_prop, v_map);
}

template <typename Mesh, typename ValueType1, typename ValueType2>
void garbage_collection(
        Mesh& _mesh,
        ExternalProperty<VH, ValueType1>& _v_prop1,
        ExternalProperty<VH, ValueType2>& _v_prop2)
{
    ISM_ASSERT(_v_prop1.size_okay(_mesh));
    ISM_ASSERT(_v_prop2.size_okay(_mesh));

    std::vector<VH> v_map;
    std::vector<HEH> h_map;
    std::vector<FH> f_map;
    garbage_collection_with_index_maps(_mesh, v_map, h_map, f_map);

    _v_prop1 = apply_index_map(_mesh, _v_prop1, v_map);
    _v_prop2 = apply_index_map(_mesh, _v_prop2, v_map);
}

}
