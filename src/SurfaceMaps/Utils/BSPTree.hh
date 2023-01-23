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
#include <SurfaceMaps/Utils/BarycentricPoint.hh>
#include <ACG_BSP/Geometry/bsp/TriangleBSPT.hh>
#include <memory>

namespace SurfaceMaps
{

struct BSPTree
{
    using OMBSP = OpenMeshTriangleBSPT<TriMesh>;
    using OMBSPptr = std::unique_ptr<OMBSP>;

    BSPTree() = default;

    BSPTree(
            const TriMesh& _mesh);

    BarycentricPoint bary(
            const Vec3d& _p,
            const TriMesh& _mesh) const;

    Vec3d project(
            const Vec3d& _p,
            const TriMesh& _mesh) const;

    SFH face(
            const Vec3d& _p,
            const TriMesh& _mesh) const;

    OMBSPptr om_bsp;
};

}
