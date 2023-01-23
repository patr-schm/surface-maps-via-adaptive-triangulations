/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt
 */

#include <SurfaceMaps/Utils/Helpers.hh>
#include <SurfaceMaps/Utils/Subdivision.hh>

namespace SurfaceMaps
{

namespace
{

BarycentricPoint mix_bary(const BarycentricPoint& _b0, const BarycentricPoint& _b1, double _t)
{
    ISM_ASSERT_EQ(_b0.heh(), _b1.heh());
    BarycentricPoint result;
    result.set_heh(_b0.heh());
    result.set_alpha_beta(
            (1.0 - _t) * _b0.alpha() + _t * _b1.alpha(),
            (1.0 - _t) * _b0.beta()  + _t * _b1.beta()
    );
    return result;
}

}

RefinementMap
subdivide_1_to_4(
        TriMesh& _mesh)
{
    return subdivide_1_to_4(_mesh, identity_refinement_map(_mesh));
}

RefinementMap
subdivide_1_to_4(
        TriMesh& _orig_mesh,
        const RefinementMap& _orig_ref_map)
{
    // Map from original edges to new vertices
    ExternalProperty<EH, VH> orig_e_to_new_v(_orig_mesh);

    // Insert original vertices
    TriMesh new_mesh;
    for (auto v : _orig_mesh.vertices())
        new_mesh.add_vertex(_orig_mesh.point(v));

    // Add new vertices
    for (auto e : _orig_mesh.edges())
        orig_e_to_new_v[e] = new_mesh.add_vertex(0.5 * _orig_mesh.point(e.v0()) + 0.5 * _orig_mesh.point(e.v1()));

    std::map<HEH, BarycentricPoint> temp_ref_map; // std::map since we're building this property as we construct the mesh

    // Add faces
    for (auto f_orig : _orig_mesh.faces())
    {
        SVH v_orig_a;
        SVH v_orig_b;
        SVH v_orig_c;
        SHEH h_orig_ab;
        SHEH h_orig_bc;
        SHEH h_orig_ca;
        handles(_orig_mesh, f_orig, v_orig_a, v_orig_b, v_orig_c, h_orig_ca, h_orig_ab, h_orig_bc);

        const SVH v_new_a = OpenMesh::make_smart(v_orig_a, new_mesh);
        const SVH v_new_b = OpenMesh::make_smart(v_orig_b, new_mesh);
        const SVH v_new_c = OpenMesh::make_smart(v_orig_c, new_mesh);

        const SVH v_new_ab = OpenMesh::make_smart(orig_e_to_new_v[h_orig_ab.edge()], new_mesh);
        const SVH v_new_bc = OpenMesh::make_smart(orig_e_to_new_v[h_orig_bc.edge()], new_mesh);
        const SVH v_new_ca = OpenMesh::make_smart(orig_e_to_new_v[h_orig_ca.edge()], new_mesh);

        const BarycentricPoint bary_a = _orig_ref_map[h_orig_ca];
        const BarycentricPoint bary_b = _orig_ref_map[h_orig_ab];
        const BarycentricPoint bary_c = _orig_ref_map[h_orig_bc];

        const BarycentricPoint bary_ab = mix_bary(bary_a, bary_b, 0.5);
        const BarycentricPoint bary_bc = mix_bary(bary_b, bary_c, 0.5);
        const BarycentricPoint bary_ca = mix_bary(bary_c, bary_a, 0.5);

        {
            new_mesh.add_face(v_new_a, v_new_ab, v_new_ca);
            temp_ref_map[new_mesh.find_halfedge(v_new_a,  v_new_ab)] = bary_ab;
            temp_ref_map[new_mesh.find_halfedge(v_new_ab, v_new_ca)] = bary_ca;
            temp_ref_map[new_mesh.find_halfedge(v_new_ca, v_new_a)]  = bary_a;
        }

        {
            new_mesh.add_face(v_new_b, v_new_bc, v_new_ab);
            temp_ref_map[new_mesh.find_halfedge(v_new_b,  v_new_bc)] = bary_bc;
            temp_ref_map[new_mesh.find_halfedge(v_new_bc, v_new_ab)] = bary_ab;
            temp_ref_map[new_mesh.find_halfedge(v_new_ab, v_new_b)]  = bary_b;
        }

        {
            new_mesh.add_face(v_new_c, v_new_ca, v_new_bc);
            temp_ref_map[new_mesh.find_halfedge(v_new_c,  v_new_ca)] = bary_ca;
            temp_ref_map[new_mesh.find_halfedge(v_new_ca, v_new_bc)] = bary_bc;
            temp_ref_map[new_mesh.find_halfedge(v_new_bc, v_new_c)]  = bary_c;
        }

        {
            new_mesh.add_face(v_new_ab, v_new_bc, v_new_ca);
            temp_ref_map[new_mesh.find_halfedge(v_new_ab, v_new_bc)] = bary_bc;
            temp_ref_map[new_mesh.find_halfedge(v_new_bc, v_new_ca)] = bary_ca;
            temp_ref_map[new_mesh.find_halfedge(v_new_ca, v_new_ab)] = bary_ab;
        }
    }

    // Construct new-to-old barycentric coordinate map
    RefinementMap ref_map(new_mesh);
    for (const auto& heh : new_mesh.halfedges())
        ref_map[heh] = temp_ref_map[heh];

    _orig_mesh = new_mesh;
    return ref_map;
}

}
