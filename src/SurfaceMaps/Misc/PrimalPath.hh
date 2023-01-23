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

struct PrimalPath
{
    std::vector<HEH> hehs;

    /// Removes successive edges that go back and forth in opposite directions.
    /// If the path is closed, _move_loop_basepoint enables removing redundant start / end elements,
    /// effectively moving the basepoint.
    template<typename Mesh>
    void shrink(const Mesh& _mesh, const bool _move_loop_basepoint = false);

    template<typename Mesh>
    bool is_closed(const Mesh& _mesh) const;

    template<typename Mesh>
    void assert_valid(const Mesh& _mesh) const;

    template<typename Mesh>
    VH vh_start(const Mesh& _mesh) const;

    template<typename Mesh>
    VH vh_end(const Mesh& _mesh) const;

    template<typename Mesh>
    double embedded_length(const Mesh& _mesh) const;

    template<typename Mesh>
    PrimalPath reversed(const Mesh& _mesh) const;
};

}
