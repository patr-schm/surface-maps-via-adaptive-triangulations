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

int count_boundaries(const TriMesh& _mesh);

bool isolated_vertices(const TriMesh& _mesh);

int euler_characteristic(const TriMesh& _mesh);

/// Genus of a closed surface.
int genus(const TriMesh& _mesh);

bool disk_topology(const TriMesh& _mesh);

bool closed_surface(const TriMesh& _mesh);

/// In constrast to genus(), this function also works for disk-topology meshes.
int count_handles(const TriMesh& _mesh);

}
