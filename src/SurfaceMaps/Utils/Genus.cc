/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt, Janis Born
 */

#include "Genus.hh"

namespace SurfaceMaps
{

int count_boundaries(const TriMesh& _mesh)
{
    int n_boundaries = 0;

    std::vector<bool> visited(_mesh.n_halfedges(), false);
    for (auto heh : _mesh.halfedges())
    {
        if (_mesh.has_halfedge_status() && _mesh.status(heh).deleted())
            continue;

        if (!_mesh.is_boundary(heh) || visited[heh.idx()])
            continue;

        ++n_boundaries;

        // Visit all halfedges of boundary loop
        for (auto hl_it = _mesh.chl_begin(heh); hl_it != _mesh.chl_end(heh); ++hl_it)
            visited[hl_it->idx()] = true;
    }

    return n_boundaries;
}

bool isolated_vertices(const TriMesh& _mesh)
{
    for (auto v : _mesh.vertices())
    {
        if (_mesh.cvoh_begin(v) == _mesh.cvoh_end(v))
            return true;
    }

    return false;
}

bool disk_topology(const TriMesh& _mesh)
{
    return !isolated_vertices(_mesh) &&
           euler_characteristic(_mesh) == 1 &&
           count_boundaries(_mesh) == 1;
}

bool closed_surface(const TriMesh& _mesh)
{
    return !isolated_vertices(_mesh) &&
           count_boundaries(_mesh) == 0;
}

int euler_characteristic(const TriMesh& _mesh)
{
    return _mesh.n_vertices() - _mesh.n_edges() + _mesh.n_faces();
}

/// Genus of a closed surface.
int genus(const TriMesh& _mesh)
{
    // Has to be a closed surface, but we removed the assertion here for speedup.
    //ISM_ASSERT(closed_surface(_mesh));
    const int euler = euler_characteristic(_mesh);
    ISM_ASSERT_EQ(euler % 2, 0);
    return 1 - euler / 2;
}

/// Other than genus(), this function also works for disk-topology meshes.
int count_handles(const TriMesh& _mesh)
{
    const int n_holes = count_boundaries(_mesh);
    const int g2 = ((int)_mesh.n_edges() - (int)_mesh.n_vertices() - (int)_mesh.n_faces() - n_holes + 2);
    ISM_ASSERT_EQ(g2 % 2, 0);

    return g2 / 2;
}

}
