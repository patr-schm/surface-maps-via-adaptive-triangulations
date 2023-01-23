/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt
 */

#include <SurfaceMaps/Utils/BSPTree.hh>
#include <TinyAD/Utils/Timer.hh>

namespace SurfaceMaps
{

BSPTree::BSPTree(
        const TriMesh& _mesh)
{
    ISM_INFO("Building BSP tree");
    TinyAD::Timer timer("Building BSP tree");

    om_bsp.reset(new OMBSP(_mesh));
    om_bsp->reserve(_mesh.n_faces());
    for (const FH fh : _mesh.faces())
        om_bsp->push_back(fh);
    om_bsp->build(10, 100); // max vertices per leaf 10, max depth 100
}


BarycentricPoint BSPTree::bary(
        const Vec3d& _p,
        const TriMesh& _mesh) const
{
    ISM_ASSERT_FINITE_MAT(_p);
    ISM_ASSERT_EQ(om_bsp->size(), _mesh.n_faces());
    const auto nn = om_bsp->nearest(to_acg3(_p));
    ISM_ASSERT(nn.handle.is_valid());
    ISM_ASSERT(_mesh.is_valid_handle(nn.handle));

    // TODO: Clip to triange!!!

    return BarycentricPoint(_p, nn.handle, _mesh);
}

Vec3d BSPTree::project(
        const Vec3d& _p,
        const TriMesh& _mesh) const
{
    return bary(_p, _mesh).point(_mesh);
}

SFH BSPTree::face(
        const Vec3d& _p,
        const TriMesh& _mesh) const
{
    return OpenMesh::make_smart(bary(_p, _mesh).fh(_mesh), _mesh);
}

}
