/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Joe Jakobi, Patrick Schmidt
 */

#include "MultiResSphereEmbedding.hh"

#include <SurfaceMaps/Utils/LocalCoordinateSystem.hh>
#include <SurfaceMaps/MultiRes/VertexOptimization.hh>
#include <SurfaceMaps/Misc/SphereEmbeddingTet.hh>
#include <SurfaceMaps/Misc/ConstantCurvatureGeometry.hh>
#include <SurfaceMaps/Misc/ConstantCurvatureEmbedding.hh>
#include <SurfaceMaps/Misc/DistortionEnergy.hh>
#include <SurfaceMaps/AdaptiveTriangulations/OptimizeSphereEmbedding.hh>
#include <TinyAD/Utils/Timer.hh>

namespace SurfaceMaps
{

namespace
{

bool is_legal_collapse(
        const TriMesh& _mesh,
        HEH _heh)
{
    if (_mesh.is_boundary(_mesh.from_vertex_handle(_heh)))
    {
        if (!_mesh.is_boundary(_mesh.edge_handle(_heh)))
            return false;
        if (_mesh.valence(_mesh.from_vertex_handle(_heh)) == 2)
            return false;
    }

    return true;

}

Vec3d slerp(
        const Vec3d& _p0,
        const Vec3d& _p1,
        const double _t)
{
    double omega = acos(_p0.dot(_p1));

    ISM_ASSERT_NEQ(_p0, _p1);

    if (omega != omega || omega < 1e-16)
    {
        ISM_ERROR("Slerp failed: omega = " << omega);
        Vec3d p_new = (1.0 - _t) * _p0 + _t * _p1;
        return project_to_model(p_new, Geometry::Spherical);
    }

    Vec3d p_new = _p0.operator*(sin(omega * (1.0 - _t)) / sin(omega)) +
            _p1.operator*(sin(_t * omega) / sin(omega));
    return project_to_model(p_new, Geometry::Spherical);
}

double compute_inscribed_circle_radius(
        const TriMesh& _mesh,
        const VH _vh_0,
        const VH _vh_1,
        const VH _vh_left,
        const VH _vh_right,
        const OpenMesh::VPropHandleT<Vec3d> _ph)
{
    std::vector<Vec3d> polygon_vertices;
    OpenMesh::SmartHalfedgeHandle heh_iter;
    for (auto heh : _mesh.voh_range(_vh_0))
    {
        if (heh.to() == _vh_left)
        {
            heh_iter = heh;
            break;
        }
    }

    ISM_ASSERT(heh_iter.is_valid());

    while (heh_iter.to() != _vh_1)
    {
        polygon_vertices.push_back(_mesh.property(_ph, heh_iter.to()));
        heh_iter = heh_iter.prev().opp();
    }

    double radius = INF_DOUBLE;
    const Vec3d p_center = _mesh.property(_ph, _vh_1);

    for (uint i = 0; i < polygon_vertices.size() - 1; ++i)
    {
        const Vec3d current_pos = polygon_vertices.at(i);
        const Vec3d next_pos = polygon_vertices.at(i + 1);

        const Vec3d n = current_pos.cross(next_pos);

//        radius = std::min(radius, asin(std::abs(p_center.dot(n))));
        radius = std::max(radius, asin(std::abs(p_center.dot(n))));
    }

    radius = std::min(radius, M_PI * 0.99);

    ISM_ASSERT_L(radius, M_PI);
    ISM_ASSERT_G(radius, 0);

    return radius;
}

bool spherical_line_search(
        const Vec3d& _p_from,
        const Vec3d& _p_to,
        const VH _vh_new,
        const MultiResSphereEmbedding& _data,
        const MultiResSphereEmbeddingSettings& _settings,
        Vec3d& _p_out)
{
    double t = 1.0;
    double best_energy = INF_DOUBLE;

    ISM_ASSERT_NOT_NAN3(_p_from);
    ISM_ASSERT_NOT_NAN3(_p_to);

    for (uint i = 0; i < 250; ++i)
    {
        Vec3d p = slerp(_p_from, _p_to, t);
        t *= 0.95;
        double current_energy = 0;

        for (auto fh : _data.mesh.vf_range(_vh_new))
        {
            VH vh_a, vh_b, vh_c;
            handles(_data.mesh, fh, vh_a, vh_b, vh_c);
            Vec3d p_a = vh_a == _vh_new ? p : _data.mesh.property(_data.ph_embedding, vh_a);
            Vec3d p_b = vh_b == _vh_new ? p : _data.mesh.property(_data.ph_embedding, vh_b);
            Vec3d p_c = vh_c == _vh_new ? p : _data.mesh.property(_data.ph_embedding, vh_c);

            ISM_ASSERT_FINITE_MAT(p);

            current_energy += sphere_distortion_energy(
                        fh, p_a, p_b, p_c,
                        _data.mesh, _data.mesh_area,
                        _settings.refinement_opt_settings.w_barrier,
                        _settings.refinement_opt_settings.w_angle,
                        _settings.refinement_opt_settings.w_area,
                        _settings.refinement_opt_settings.max_edge_length_degrees);

        }

        if (current_energy < best_energy)
        {
            best_energy = current_energy;
            _p_out = p;
        }
    }

    ISM_EXPECT_FINITE(best_energy);
    if (!isfinite(best_energy))
        return false;

    return true;
}

bool compute_position_in_kernel(
        MultiResSphereEmbedding& _data,
        const VH _vh_0,
        const VH _vh_1,
        const VH _vh_left,
        const VH _vh_right,
        const MultiResSphereEmbeddingSettings& _settings,
        bool _try_other_vh_1,
        Vec3d& _p_in_kernel)
{
    const Vec3d p_l = _data.mesh.property(_data.ph_embedding, _vh_left);
    const Vec3d p_r = _data.mesh.property(_data.ph_embedding, _vh_right);
    const Vec3d p_1 = _data.mesh.property(_data.ph_embedding, _vh_1);

    // Check if angle between (v_l, v_1) and (v_1, v_r) is greater than 180°
    // For this, calc normal n of (v_1, v_l). Angle is >=180° iff n^T v_r <= 0
    const Vec3d n_l = p_1.cross(p_l).normalized();
    const Vec3d n_r = p_r.cross(p_1).normalized();

    Vec3d dir = n_l + n_r;
    if (dir.norm() < 1e-3)
    {
        ISM_DEBUG_OUT("Use midpoint to determine dir");
        const bool greater_pi = (n_l.dot(p_r) <= 0);

        const double dist_l = distance(p_1, p_l, Spherical);
        const double dist_r = distance(p_1, p_r, Spherical);

        ISM_EXPECT_G(dist_l, 0);
        ISM_EXPECT_G(dist_r, 0);

        // Compute weighted midpoint between p_l and p_r
        Vec3d p_mid;
        if (dist_l < dist_r)
        {
            const Vec3d p = slerp(p_1, p_r, dist_l/dist_r);
            p_mid = slerp(p_l, p, 0.5);
        }
        else
        {
            const Vec3d p = slerp(p_1, p_l, dist_r/dist_l);
            p_mid = slerp(p_r, p, 0.5);
        }

        dir = greater_pi ? p_1 - p_mid  : p_mid - p_1;
    }
    dir.normalize();
    ISM_ASSERT_NOT_NAN3(dir);

    const double radius_inscribed_circle = compute_inscribed_circle_radius(_data.mesh, _vh_0, _vh_1, _vh_left, _vh_right, _data.ph_embedding);

    Vec3d p_to = p_1 + dir;
    p_to =  project_to_model(p_to, Spherical);

    const Vec3d p_target = slerp(p_1, p_to,  radius_inscribed_circle / distance(p_1, p_to, Spherical));

    if (!spherical_line_search(p_1, p_target, _vh_0, _data, _settings, _p_in_kernel))
    {
        if (!_try_other_vh_1)
        {
            ISM_ERROR("ERROR in compute_position_in_kernel");
            return false;
        }
        for (SHEH heh : _data.mesh.voh_range(_vh_0))
        {
            return compute_position_in_kernel(_data, _vh_0, heh.to(), heh.next().to(),
                                              heh.opp().next().to(), _settings, false, _p_in_kernel);
        }

        ISM_ERROR("No position found!");
        return false;
    }

    for (auto heh : _data.mesh.voh_range(_vh_0))
    {
        ISM_ASSERT(!flipped_or_degenerate(
                        _p_in_kernel,
                        _data.mesh.property(_data.ph_embedding, heh.to()),
                        _data.mesh.property(_data.ph_embedding, heh.next().to()),
                        Geometry::Spherical));
    }

    return true;
}

bool refine_independent_set(
        MultiResSphereEmbedding& _data,
        const ProgressiveMesh& _pm,
        const MultiResSphereEmbeddingSettings& _settings)
{
    // Refine Mesh
    std::vector<VH> vhs_optimize;
    while (vhs_optimize.size() < _data.independent_set_sizes.at(_data.current_independent_set))
    {
        // Insert vertex (connectivity)
        undo_decimation_step(_data.mesh, *_data.rec_it);
        _data.mesh_area = total_area(_data.mesh); // TODO: Why this? Is looks slow

        // Compute initial vertex position
        // and write to property
        if (!compute_position_in_kernel(
                    _data, _data.rec_it->vh0, _data.rec_it->vh1, _data.rec_it->vhl, _data.rec_it->vhr,
                    _settings, true, _data.mesh.property(_data.ph_embedding, _data.rec_it->vh0)))
            return false;

        // Sanity check
        for(auto heh : _data.mesh.voh_range(_data.rec_it->vh0))
        {
            ISM_ASSERT(!flipped_or_degenerate(_data.mesh.property(_data.ph_embedding, heh.from()),
                                              _data.mesh.property(_data.ph_embedding, heh.to()),
                                              _data.mesh.property(_data.ph_embedding, heh.next().to()),
                                              Spherical));
        }

        vhs_optimize.push_back(_data.rec_it++->vh0);
    }
    --_data.current_independent_set;

    // Optimize newly inserted vertices
    if (_settings.parallel_optimization)
    {
        #pragma omp parallel for schedule(dynamic)
        for (uint i=  0; i < vhs_optimize.size(); ++i)
        {
            optimize_vertex(_data, vhs_optimize.at(i),
                            _settings.initialization_max_iters, _settings);
        }
    }
    else
    {
        for (uint i=  0; i < vhs_optimize.size(); ++i)
            optimize_vertex(_data, vhs_optimize.at(i),
                            _settings.initialization_max_iters, _settings);
    }

    ISM_INFO("Inserted " << _data.mesh.n_vertices() << " of " << _pm.orig_to_prog_idx.size() << " vertices.");
    return true;
}

void global_optimization(
        MultiResSphereEmbedding& _data,
        const MultiResSphereEmbeddingSettings& _settings,
        std::function<void(const TriMesh&, const ExternalProperty<VH, Vec3d>&)> _callback)
{
    // Determine n_rounds
    int n_rounds;
    if (_settings.use_phases)
    {
        if ((int)_data.mesh.n_vertices() < _settings.phase_one_vertex_thresh)
            n_rounds = _settings.global_iters_phase_one;
        else
            n_rounds = _settings.global_iters_phase_two;
    }
    else
    {
        n_rounds =  std::pow(1.0 - 1.0 * _data.mesh.n_vertices() / _data.n_vertices_total, _settings.iters_exponent)
                * _settings.iters_exponentional;

        n_rounds = std::max(n_rounds, _settings.iters_minimal);
    }

    ISM_DEBUG_OUT("Starting " << n_rounds << " iterations of global optimization");
    TinyAD::Timer timer("Global sphere optimization");

    // Global Optimization
    ExternalProperty<VH, Vec3d> embedding(_data.mesh);
    std::vector<VH> constraints;

    #pragma omp parallel for
    for (uint i = 0; i < _data.mesh.n_vertices(); ++i)
        embedding[VH(i)] = _data.mesh.property(_data.ph_embedding, VH(i));

    optimize_sphere_embedding(_data.mesh, embedding, n_rounds, _settings.refinement_opt_settings, [&] ()
    {
        _callback(_data.mesh, embedding);
    });

    #pragma omp parallel for
    for (uint i = 0; i < _data.mesh.n_vertices(); ++i)
        _data.mesh.property(_data.ph_embedding, VH(i)) = embedding[VH(i)];
}

bool refine_and_optimize(
        MultiResSphereEmbedding& _data,
        const ProgressiveMesh& _pm,
        const MultiResSphereEmbeddingSettings& _settings,
        std::function<void(const TriMesh&, const ExternalProperty<VH, Vec3d>&)> _callback)
{
    TinyAD::Timer timer("Embedding refinement and optimization");

    while ((int)_data.mesh.n_vertices() != _data.n_vertices_total)
    {
        if (!refine_independent_set(_data, _pm, _settings))
            return false;

        global_optimization(_data, _settings, _callback);
    }

    {
        ISM_DEBUG_OUT("Starting " << _settings.global_iters_final << " final iterations of global optimization");
        TinyAD::Timer timer_final("Global sphere optimization");

        // Global Optimization
        ExternalProperty<VH, Vec3d> embedding(_data.mesh);

        #pragma omp parallel for
        for (uint i = 0; i < _data.mesh.n_vertices(); ++i)
            embedding[VH(i)] = _data.mesh.property(_data.ph_embedding, VH(i));

        optimize_sphere_embedding(_data.mesh, embedding, _settings.global_iters_final, _settings.refinement_opt_settings, [&] ()
        {
            _callback(_data.mesh, embedding);
        });

        #pragma omp parallel for
        for (uint i = 0; i < _data.mesh.n_vertices(); ++i)
            _data.mesh.property(_data.ph_embedding, VH(i)) = embedding[VH(i)];
    }

    return true;
}

bool multi_res_sphere_embedding(
        const TriMesh& _mesh,
        const MultiResSphereEmbeddingSettings& _settings,
        ExternalProperty<VH, Vec3d>& _embedding,
        std::function<void(const TriMesh&, const ExternalProperty<VH, Vec3d>&)> _callback)
{
    TinyAD::Timer timer("Multi Res Parametrization");

    // Decimate genus 0 mesh to tetrahedron
    TriMesh mesh_tmp = _mesh;
    std::vector<uint> independent_set_sizes;
    const uint target_vertices = 4;
    ProgressiveMesh progressive_mesh = compute_progressive_mesh(
                mesh_tmp,
                target_vertices,
                _settings.decimation_settings,
                independent_set_sizes,
                [&] (HEH _heh) { return is_legal_collapse(_mesh, _heh); });
    ISM_ASSERT_EQ(mesh_tmp.n_vertices(), target_vertices);

    // Embed tetrahedron on sphere
    ExternalProperty<VH, Vec3d> embedding_ext_prop;
    sphere_embedding_tet(mesh_tmp, EH(0), embedding_ext_prop);
    OpenMesh::VPropHandleT<Vec3d> ph_embedding;
    mesh_tmp.add_property(ph_embedding);
    for (VH vh : mesh_tmp.vertices())
        mesh_tmp.property(ph_embedding, vh) = embedding_ext_prop[vh];

    _callback(mesh_tmp, embedding_ext_prop);

    // Refine and optimize mesh on sphere
    MultiResSphereEmbedding data(mesh_tmp, _mesh.n_vertices(), ph_embedding, progressive_mesh, independent_set_sizes);
    if (!refine_and_optimize(data, progressive_mesh, _settings, _callback))
        return false;

    _embedding.init(_mesh);
    for (auto v : _mesh.vertices())
        _embedding[v] = data.mesh.property(data.ph_embedding, VH(data.pm.orig_to_prog_idx.at(v.idx())));

    return true;
}

}

ExternalProperty<VH, Vec3d> multi_res_sphere_embedding(
        const TriMesh& _mesh,
        const MultiResSphereEmbeddingSettings& _settings,
        std::function<void(const TriMesh&, const ExternalProperty<VH, Vec3d>&)> _callback)
{
    ExternalProperty<VH, Vec3d> embedding;
    bool success = multi_res_sphere_embedding(_mesh, _settings, embedding, _callback);

    if (success)
        return embedding;

    if (_settings.try_hard)
    {
        // Try a bunch of settings
        MultiResSphereEmbeddingSettings try_hard_settings = _settings;

        try_hard_settings.refinement_opt_settings.w_barrier = 1000.0;
        try_hard_settings.refinement_opt_settings.w_angle = 100.0;
        try_hard_settings.refinement_opt_settings.w_area = 1.0;
        try_hard_settings.initialization_max_iters = 20;
        try_hard_settings.global_iters_phase_one = 50;
        try_hard_settings.global_iters_phase_two = 50;

        ISM_WARNING("Multi res sphere embedding failed. Trying again with different settings.");
        success = multi_res_sphere_embedding(_mesh, try_hard_settings, embedding, _callback);
        if (success)
            return embedding;

        try_hard_settings.refinement_opt_settings.w_angle = 0.0;
        try_hard_settings.refinement_opt_settings.w_area = 0.0;

        ISM_WARNING("Multi res sphere embedding failed. Trying again with different settings.");
        success = multi_res_sphere_embedding(_mesh, try_hard_settings, embedding, _callback);
        if (success)
            return embedding;
    }

    ISM_ERROR_throw("Multi res sphere embedding failed.");
}

}
