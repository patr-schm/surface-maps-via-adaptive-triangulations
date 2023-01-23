/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt, Janis Born
 */

#include "PrimalPath.hh"

namespace SurfaceMaps
{

template<typename Mesh>
void PrimalPath::shrink(
        const Mesh& _mesh,
        const bool _move_loop_basepoint)
{
    bool good = false;
    do
    {
        // (Optional) cancellations that affect loop start / endpoints.
        if (_move_loop_basepoint && is_closed(_mesh))
        {
            const auto& heh_start = hehs.front();
            const auto& heh_end_opp = _mesh.opposite_halfedge_handle(hehs.back());
            if (heh_start == heh_end_opp)
            {
                hehs.erase(hehs.begin());
                hehs.pop_back();
                good = false;
                ISM_ASSERT(is_closed(_mesh));
                continue;
            }
        }
        // Interior cancellations
        for (int i0 = 0; i0 < hehs.size() - 1; ++i0)
        {
            const int i1 = i0 + 1;
            const auto& heh0 = hehs[i0];
            const auto& heh1 = hehs[i1];
            const auto& heh0_opp = _mesh.opposite_halfedge_handle(heh0);
            if (heh1 == heh0_opp)
            {
                hehs.erase(hehs.begin() + i0, hehs.begin() + i1 + 1);
                good = false;
                break;
            }
        }
        good = true;
    }
    while (!good);
    assert_valid(_mesh);
}
template void PrimalPath::shrink<TriMesh>(const TriMesh&, const bool);
template void PrimalPath::shrink<PolyMesh>(const PolyMesh&, const bool);

template<typename Mesh>
bool PrimalPath::is_closed(const Mesh& _mesh) const
{
    if (hehs.size() == 0)
        return true;
    else
        return vh_start(_mesh) == vh_end(_mesh);
}
template bool PrimalPath::is_closed<TriMesh>(const TriMesh&) const;
template bool PrimalPath::is_closed<PolyMesh>(const PolyMesh&) const;

template<typename Mesh>
void PrimalPath::assert_valid(const Mesh& _mesh) const
{
    if (hehs.size() < 2)
        return;

    for (int i0 = 0; i0 < hehs.size() - 1; ++i0)
    {
        const int i1 = (i0 + 1) % hehs.size();
        const auto& heh0 = hehs[i0];
        const auto& heh1 = hehs[i1];
        ISM_ASSERT(heh0.is_valid());
        ISM_ASSERT(heh1.is_valid());
        const auto& heh0_to   = _mesh.to_vertex_handle(heh0);
        const auto& heh1_from = _mesh.from_vertex_handle(heh1);
        ISM_ASSERT(heh0_to == heh1_from);
    }
}
template void PrimalPath::assert_valid<TriMesh>(const TriMesh&) const;
template void PrimalPath::assert_valid<PolyMesh>(const PolyMesh&) const;

template<typename Mesh>
VH PrimalPath::vh_start(const Mesh& _mesh) const
{
    ISM_ASSERT(!hehs.empty());
    return _mesh.from_vertex_handle(hehs.front());
}
template VH PrimalPath::vh_start<TriMesh>(const TriMesh&) const;
template VH PrimalPath::vh_start<PolyMesh>(const PolyMesh&) const;

template<typename Mesh>
VH PrimalPath::vh_end(const Mesh& _mesh) const
{
    ISM_ASSERT(!hehs.empty());
    return _mesh.to_vertex_handle(hehs.back());
}
template VH PrimalPath::vh_end<TriMesh>(const TriMesh&) const;
template VH PrimalPath::vh_end<PolyMesh>(const PolyMesh&) const;

template<typename Mesh>
double PrimalPath::embedded_length(const Mesh& _mesh) const
{
    double result = 0.0;
    for (const auto& heh : hehs)
        result += _mesh.calc_edge_length(heh);
    return result;
}
template double PrimalPath::embedded_length<TriMesh>(const TriMesh&) const;
template double PrimalPath::embedded_length<PolyMesh>(const PolyMesh&) const;

template<typename Mesh>
PrimalPath PrimalPath::reversed(const Mesh& _mesh) const
{
    PrimalPath result = *this;
    std::reverse(begin(result.hehs), end(result.hehs));
    for (auto& heh : result.hehs)
    {
        heh = _mesh.opposite_halfedge_handle(heh);
    }
    return result;
}
template PrimalPath PrimalPath::reversed<TriMesh>(const TriMesh&) const;
template PrimalPath PrimalPath::reversed<PolyMesh>(const PolyMesh&) const;

}
