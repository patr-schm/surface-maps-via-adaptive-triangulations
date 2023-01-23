/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt, Janis Born
 */

#include "Dijkstra.hh"

#include <SurfaceMaps/Utils/Helpers.hh>
#include <queue>

namespace SurfaceMaps
{

PrimalPath find_primal_path(
        const VH _vh_source,
        const VH _vh_target,
        const TriMesh& _mesh,
        const ExternalProperty<EH, bool>& _blocked_edges,
        ExternalProperty<VH, bool> _blocked_vertices)
{
    PrimalPath path;

    if (_vh_source == _vh_target)
        return path;

    struct DijkstraEdge
    {
        DijkstraEdge(const HEH _heh, const double _distance)
            : heh(_heh), distance(_distance) { }

        // Higher priority if lower distance
        bool operator<(const DijkstraEdge& other) const { return distance > other.distance; }

        HEH heh;
        double distance; // from-vertex to to-vertex
    };

    // Block all vertices incident to blocked edges
    for (auto e : _mesh.edges())
    {
        if (_blocked_edges[e])
        {
            _blocked_vertices[e.v0()] = true;
            _blocked_vertices[e.v1()] = true;
        }
    }

    ExternalProperty<VH, bool> visited(_mesh, false);
    ExternalProperty<VH, HEH> pred(_mesh, HEH(-1)); // heh pointing to vertex
    std::priority_queue<DijkstraEdge> queue;

    auto enqueue = [&] (const HEH heh, const double dist_curr)
    {
        const VH vh_to = _mesh.to_vertex_handle(heh);
        if (visited[vh_to])
            return;

        if (_blocked_vertices[vh_to] && vh_to != _vh_target)
            return;

        const double length = _mesh.calc_edge_length(heh);
        queue.push(DijkstraEdge(heh, dist_curr + length));
    };

    for (const HEH heh : _mesh.voh_range(_vh_source))
        enqueue(heh, 0.0);

    // Dijkstra
    VH vh_curr(-1);
    while (!queue.empty())
    {
        const DijkstraEdge to_curr = queue.top();
        queue.pop();

        // Already found a shorter path?
        vh_curr = _mesh.to_vertex_handle(to_curr.heh);
        if (visited[vh_curr])
            continue;

        visited[vh_curr] = true;
        pred[vh_curr] = to_curr.heh;

        // Reached target?
        if (vh_curr == _vh_target)
            break;

        // Add neighbors
        for (const HEH heh : _mesh.voh_range(vh_curr))
            enqueue(heh, to_curr.distance);
    }

    if (vh_curr != _vh_target)
    {
        ISM_WARNING("Primal Dijkstra did not find target vertex.");
        return path;
    }

    // Trace path from target to source.
    while (vh_curr != _vh_source)
    {
        const HEH heh_pred = pred[vh_curr];
        ISM_ASSERT(_mesh.is_valid_handle(heh_pred));
        path.hehs.push_back(heh_pred);

        vh_curr = _mesh.from_vertex_handle(heh_pred);
    }

    std::reverse(path.hehs.begin(), path.hehs.end());
    return path;
}

}
