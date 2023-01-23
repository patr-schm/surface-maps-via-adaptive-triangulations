/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt, Janis Born
 */
#pragma once

#include <SurfaceMaps/Misc/PrimalPath.hh>
#include <SurfaceMaps/Types.hh>

namespace SurfaceMaps
{

PrimalPath find_primal_path(
        const VH _vh_source,
        const VH _vh_target,
        const TriMesh& _mesh,
        const ExternalProperty<EH, bool>& _blocked_edges,
        ExternalProperty<VH, bool> _blocked_vertices);

inline PrimalPath find_primal_path(
        const VH _vh_source,
        const VH _vh_target,
        const TriMesh& _mesh)
{
    const ExternalProperty<EH, bool> _blocked_edges(_mesh, false);
    const ExternalProperty<VH, bool> _blocked_vertices(_mesh, false);
    return find_primal_path(_vh_source, _vh_target, _mesh, _blocked_edges, _blocked_vertices);
}

}
