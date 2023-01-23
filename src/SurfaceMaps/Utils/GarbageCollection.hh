/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt
 */
#pragma once

#include <SurfaceMaps/Utils/ExternalProperty.hh>

namespace SurfaceMaps
{

/**
 * Garbage-collect mesh and keep vertex property valid.
 */
template <typename Mesh, typename ValueType>
void garbage_collection(
        Mesh& _mesh,
        ExternalProperty<VH, ValueType>& _v_prop);

/**
 * Garbage-collect mesh and keep two vertex properties valid.
 */
template <typename Mesh, typename ValueType1, typename ValueType2>
void garbage_collection(
        Mesh& _mesh,
        ExternalProperty<VH, ValueType1>& _v_prop1,
        ExternalProperty<VH, ValueType2>& _v_prop2);

}

#include <SurfaceMaps/Utils/GarbageCollectionImpl.hh>
