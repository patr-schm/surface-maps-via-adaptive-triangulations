/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Janis Born, Joe Jakobi, Patrick Schmidt
 */

#include "ProgressiveMesh.hh"

#include <TinyAD/Utils/Timer.hh>
#include <queue>

namespace SurfaceMaps
{

namespace
{

double min_area_priority(const TriMesh& _mesh, const HEH _heh)
{
    const VH vh_to  = _mesh.to_vertex_handle(_heh);
    const Vec3d p_A = _mesh.point(vh_to);
    double min_area = INF_DOUBLE;
    for (auto heh : _mesh.voh_range(_mesh.from_vertex_handle(_heh)))
    {
        if  (heh.to() == vh_to || heh.next().to() == vh_to)
            continue;

        const Vec3d p_B = _mesh.point(heh.to());
        const Vec3d p_C = _mesh.point(heh.next().to());
        const Vec3d v_ab = p_B - p_A;
        const Vec3d v_ac = p_C - p_A;
        const double area = std::max(v_ab.cross(v_ac).norm() / 2, 1e-16);
        min_area = std::min(min_area, area);
    }
    return 1.0 / (min_area * min_area);
}

double random_priority(const HEH _heh)
{
    return _heh.idx();
}

double dihedral_angle_change_priority(
        const TriMesh& _mesh,
        const HEH _heh)
{
    // Priority: angle change of face normals
    const auto vh_from = _mesh.from_vertex_handle(_heh);
    const auto vh_to = _mesh.to_vertex_handle(_heh);
    const auto& p_from = _mesh.point(vh_from);
    const auto& p_to = _mesh.point(vh_to);
    double max_angle = 0.0;
    for (const auto& oheh : _mesh.voh_range(vh_from))
    {
        if (_mesh.is_boundary(oheh))
        {
            continue;
        }
        const auto vh1 = _mesh.to_vertex_handle(oheh);
        const auto vh2 = _mesh.opposite_vh(oheh);
        if (vh1 == vh_to || vh2 == vh_to)
        {
            continue;
        }
        const auto& p1 = _mesh.point(vh1);
        const auto& p2 = _mesh.point(vh2);
        const auto n_before = ((p1 - p_from).cross(p2 - p_from)).normalized();
        const auto n_after = ((p1 - p_to).cross(p2 - p_to)).normalized();
        const auto cos_angle = n_before.dot(n_after);
        const auto angle = std::acos(cos_angle);
        max_angle = std::max(max_angle, angle);
    }
    return max_angle;
}

double uniform_priority(
        const TriMesh& _mesh,
        const HEH _heh)
{
    return _mesh.calc_edge_length(_heh);
}


MatXd get_quadric(
        const TriMesh& _mesh,
        const VH _vh)
{
    Vec3d v = _mesh.point(_vh);
    MatXd Q = MatXd::Zero(4,4);

    for (FH fh : _mesh.vf_range(_vh))
    {
        Vec3d n = _mesh.calc_face_normal(fh);
        double a = n.x();
        double b = n.y();
        double c = n.z();

        double d = -n.dot(v);
        MatXd Q_i(4,4);
        Q_i << a*a, a*b, a*c, a*d,
             a*b, b*b, b*c, b*d,
             a*c, b*c, c*c, c*d,
             a*d, b*d, c*d, d*d;
        Q += Q_i;
    }

    return Q;
}

double error_quadrics_priority(
        const TriMesh& _mesh,
        const HEH _heh,
        OpenMesh::VPropHandleT<MatXd>& _ph_quadrics)
{
    const VH vh1 = _mesh.from_vertex_handle(_heh);
    const VH vh2 = _mesh.to_vertex_handle(_heh);

    MatXd Q1 = _mesh.property(_ph_quadrics, vh1);
    MatXd Q2 = _mesh.property(_ph_quadrics, vh2);

    TriMesh::Point p = _mesh.point(vh2);
    Vec4d v(p.x(), p.y(), p.z(), 1);

    double res = v.transpose() * (Q1 + Q2) * v;

//    ISM_EXPECT_GEQ(res, 0);
    ISM_EXPECT_GEQ(res, -1e-12);

    return std::max(res, 0.0);
}

double get_angle(
        const Vec3d& _e0,
        const Vec3d& _e1)
{
    const double dot = _e0.dot(_e1);
    return acos(std::clamp(dot / (_e0.norm() * _e1.norm()), -1.0, 1.0));
}

double angle_defect(
        const TriMesh& _mesh,
        const VH _vh)
{
    double total_angle = 0;
    for (HEH heh: _mesh.vih_range(_vh))
    {
        if (_mesh.is_boundary(heh))
            continue;

        total_angle += _mesh.calc_sector_angle(heh);
    }

    if (_mesh.is_boundary(_vh))
        return M_PI - total_angle;
    else
        return 2.0 * M_PI - total_angle;
}

double angle_defect_to_vertex(
        const TriMesh& _mesh,
        const VH _vh_to,
        const VH _vh_center,
        const VH _vh_left,
        const VH _vh_right)
{
    const TriMesh::Point p0 = _mesh.point(_vh_to);
    TriMesh::Point p1;
    TriMesh::Point p2;

    double total_angle = 0.0;

    if (_mesh.is_valid_handle(_vh_right))
    {
        // sum up angle in the section where mesh connectivity is changed
        const auto heh_right = _mesh.find_halfedge(_vh_center, _vh_right);
        auto heh_iter = heh_right;
        while (heh_iter.to() != _vh_left && !heh_iter.opp().is_boundary())
        {
            ISM_ASSERT_EQ(heh_iter.from(), _vh_center)
            p1 = _mesh.point(heh_iter.to());
            heh_iter = heh_iter.opp();
            heh_iter = heh_iter.next();
            ISM_ASSERT_EQ(heh_iter.from(), _vh_center)
            p2 = _mesh.point(heh_iter.to());

            const Vec3d e0 = p1 - p0;
            const Vec3d e1 = p2 - p0;

            if (p1 == p0 || p2 == p0)
                continue;

            total_angle += get_angle(e0, e1);
        }

        ISM_ASSERT(!std::isnan(total_angle));

        // sum up angle of the unchanged part
        heh_iter = heh_right;
        heh_iter = heh_iter.next();
        heh_iter = heh_iter.opp();
        while (heh_iter.to() != _vh_left && !heh_iter.is_boundary())
        {
            p1 = _mesh.point(heh_iter.to());
            heh_iter = heh_iter.prev();
            heh_iter = heh_iter.opp();
            p2 = _mesh.point(heh_iter.to());

            const Vec3d e0 = p1 - p0;
            const Vec3d e1 = p2 - p0;

            if (p0 == p1 || p0 == p2)
                continue;

            total_angle += get_angle(e0, e1);
        }

        if (heh_iter.to() != _vh_left && _vh_left.is_valid())
        {
            heh_iter = _mesh.find_halfedge(_vh_to, _vh_left);
            while (!heh_iter.opp().is_boundary())
            {
                p1 = _mesh.point(heh_iter.to());
                ISM_ASSERT_EQ(heh_iter.from(), _vh_to);
                heh_iter = heh_iter.opp();
                ISM_ASSERT(!heh_iter.is_boundary());
                heh_iter = heh_iter.next();
                ISM_ASSERT_EQ(heh_iter.from(), _vh_to);
                p2 = _mesh.point(heh_iter.to());

                const Vec3d e0 = p1 - p0;
                const Vec3d e1 = p2 - p0;

                if (p0 == p1 || p0 == p2)
                    continue;

                total_angle += get_angle(e0, e1);
            }
        }
    }
    else if (_mesh.is_valid_handle(_vh_left)) // _vh_right is unvalid
    {
        // sum up angle in the section where mesh connectivity is changed
        const auto heh_left = _mesh.find_halfedge(_vh_center, _vh_left);
        auto heh_iter = heh_left;
        while(!heh_iter.is_boundary())
        {
            p1 = _mesh.point(heh_iter.to());
            heh_iter = heh_iter.prev();
            heh_iter = heh_iter.opp();
            p2 = _mesh.point(heh_iter.to());

            const Vec3d e0 = p1 - p0;
            const Vec3d e1 = p2 - p0;

            if (p1 == p0 || p2 == p0)
                continue;

            total_angle += get_angle(e0, e1);
        }

        // sum up angle of the unchanged part
        heh_iter = heh_left.opp();
        heh_iter = heh_iter.prev();
        while (!heh_iter.opp().is_boundary())
        {
            p1 = _mesh.point(heh_iter.to());
            heh_iter = heh_iter.opp();
            heh_iter = heh_iter.next();
            p2 = _mesh.point(heh_iter.to());

            const Vec3d e0 = p1 - p0;
            const Vec3d e1 = p2 - p0;

            if (p0 == p1 || p0 == p2)
                continue;

            total_angle += get_angle(e0, e1);
        }
    }
    else
    {
        ISM_ERROR_throw("left and right vertices were not valid");
    }

    ISM_ASSERT_NOT_NAN(total_angle);

    if (_mesh.is_boundary(_vh_to))
        return M_PI - total_angle;
    else
        return 2.0 * M_PI - total_angle;
}

double angle_defect_left_vertex(
        const TriMesh& _mesh,
        const VH _vh_left,
        const VH _vh_center,
        const VH _vh_to)
{
    double total_angle = 0;

    const TriMesh::Point p0 = _mesh.point(_vh_left);
    for (auto heh: _mesh.voh_range(_vh_left))
    {
        if (heh.is_boundary() || heh.to() == _vh_center)
            continue;

        TriMesh::Point p1;
        TriMesh::Point p2;
        if (heh.next().to() == _vh_center)
        {
            p1 = _mesh.point(heh.to());
            p2 = _mesh.point(_vh_to);
        }
        else
        {
            p1 = _mesh.point(heh.to());
            p2 = _mesh.point(heh.next().to());
        }

        const Vec3d e0 = p1 - p0;
        const Vec3d e1 = p2 - p0;

        if (p1 == p0 || p2 == p0)
            continue;

        total_angle += get_angle(e0, e1);
    }
    ISM_ASSERT_NOT_NAN(total_angle);
    if (_mesh.is_boundary(_vh_left))
        return M_PI - total_angle;
    else
        return 2.0 * M_PI - total_angle;
}

double angle_defect_right_vertex(
        const TriMesh& _mesh,
        const VH _vh_right,
        const VH _vh_center,
        const VH _vh_to)
{
    double total_angle = 0;

    const TriMesh::Point p0 = _mesh.point(_vh_right);
    for (auto heh: _mesh.voh_range(_vh_right))
    {
        if (heh.is_boundary() || heh.to() == _vh_center)
            continue;

        TriMesh::Point p1;
        TriMesh::Point p2;
        if (heh.to() == _vh_to)
        {
            auto heh_iter = heh.next();
            heh_iter = heh_iter.next();
            heh_iter = heh_iter.opp();
            heh_iter = heh_iter.next();

            p1 = _mesh.point(heh.to());
            p2 = _mesh.point(heh_iter.to());
        }
        else
        {
            p1 = _mesh.point(heh.to());
            p2 = _mesh.point(heh.next().to());
        }

        const Vec3d e0 = p1 - p0;
        const Vec3d e1 = p2 - p0;

        if (p1 == p0 || p2 == p0)
            continue;

        total_angle += get_angle(e0, e1);
    }

    ISM_ASSERT_NOT_NAN(total_angle);

    if (_mesh.is_boundary(_vh_right))
        return M_PI - total_angle;
    else
        return 2.0 * M_PI - total_angle;
}

double angle_defect_other_vertex(
        const TriMesh& _mesh,
        const VH _vh,
        const VH _vh_center,
        const VH _vh_to)
{
    const TriMesh::Point p0 = _mesh.point(_vh);

    double total_angle = 0;

    for (auto heh : _mesh.voh_range(_vh))
    {
        if (heh.is_boundary())
            continue;

        TriMesh::Point p1;
        if (heh.to() == _vh_center)
            p1 = _mesh.point(_vh_to);
        else
            p1 = _mesh.point(heh.to());

        TriMesh::Point p2;

        if (heh.next().to() == _vh_center)
            p2 = _mesh.point(_vh_to);
        else
            p2 = _mesh.point(heh.next().to());

        const Vec3d e0 = p1 - p0;
        const Vec3d e1 = p2 - p0;

        if (p1 == p0 || p2 == p0)
            continue;

        total_angle += get_angle(e0, e1);
    }

    ISM_ASSERT_NOT_NAN(total_angle);

    if (_mesh.is_boundary(_vh))
        return M_PI - total_angle;
    else
        return 2.0 * M_PI - total_angle;
}

double angle_defect_change_without_memory_priority(
        TriMesh& _mesh,
        const HEH _heh)
{
    const VH vh_center = _mesh.from_vertex_handle(_heh);
    const double center_ad_change = std::abs(angle_defect(_mesh, vh_center));

    const VH vh_to = _mesh.to_vertex_handle(_heh);
    VH vh_left = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(_heh));
    VH vh_right = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(_heh)));

    if (_mesh.is_boundary(_heh))
        vh_left = VH(-1);
    if (_mesh.is_boundary(_mesh.opposite_halfedge_handle(_heh)))
        vh_right = VH(-1);

    double one_ring_ad_change = 0;

    for (VH vh: _mesh.vv_range(vh_center))
    {
        const double ad_before_collapse = angle_defect(_mesh, vh);
        double ad_after_collapse = 0;

        if (vh == vh_to)
            ad_after_collapse = angle_defect_to_vertex(_mesh, vh_to, vh_center, vh_left, vh_right);
        else if (vh == vh_left)
            ad_after_collapse = angle_defect_left_vertex(_mesh, vh_left, vh_center, vh_to);
        else if (vh == vh_right)
            ad_after_collapse = angle_defect_right_vertex(_mesh, vh_right, vh_center, vh_to);
        else
            ad_after_collapse = angle_defect_other_vertex(_mesh, vh, vh_center, vh_to);

        one_ring_ad_change += std::abs(ad_before_collapse - ad_after_collapse);
    }

    return center_ad_change + one_ring_ad_change;
}

double angle_defect_change_with_memory_priority(
        TriMesh& _mesh,
        const HEH _heh)
{
    // Add Property to store initial angle_defect
    OpenMesh::VPropHandleT<double> ph_initial_angle_defect;
    //#pragma omp critical
    if (!_mesh.get_property_handle(ph_initial_angle_defect, "Initial Angle Defect Property Handle"))
    {
        _mesh.add_property(ph_initial_angle_defect, "Initial Angle Defect Property Handle");

        // init properties
        for (VH vh : _mesh.vertices())
            _mesh.property(ph_initial_angle_defect, vh) = angle_defect(_mesh, vh);
    }

    // Define some vertices
    const VH vh_center = _mesh.from_vertex_handle(_heh);

    // Store angle defect of center vertex
    const double center_ad_change = std::abs(_mesh.property(ph_initial_angle_defect, vh_center));

    const VH vh_to = _mesh.to_vertex_handle(_heh);

    VH vh_left = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(_heh));
    VH vh_right = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(_mesh.opposite_halfedge_handle(_heh)));

    if (_mesh.is_boundary(_heh))
        vh_left = VH(-1);
    if (_mesh.is_boundary(_mesh.opposite_halfedge_handle(_heh)))
        vh_right = VH(-1);

    double one_ring_ad_change = 0;

    for (VH vh: _mesh.vv_range(vh_center))
    {
        double ad_before_collapse = _mesh.property(ph_initial_angle_defect, vh);

        double ad_after_collapse = 0;

        if (vh == vh_to)
            ad_after_collapse = angle_defect_to_vertex(_mesh, vh_to, vh_center, vh_left, vh_right);
        else if (vh == vh_left)
            ad_after_collapse = angle_defect_left_vertex(_mesh, vh_left, vh_center, vh_to);
        else if (vh == vh_right)
            ad_after_collapse = angle_defect_right_vertex(_mesh, vh_right, vh_center, vh_to);
        else
            ad_after_collapse = angle_defect_other_vertex(_mesh, vh, vh_center, vh_to);

        one_ring_ad_change += std::abs(ad_after_collapse - ad_before_collapse);
    }

    double result = center_ad_change + one_ring_ad_change;

    return result;
}

bool area_zero_after_collapse(
        TriMesh& _mesh,
        HEH _heh,
        const ProgressiveMeshOptions& _options)
{
    const VH vh_A = _mesh.to_vertex_handle(_heh);
    auto  p_A = _mesh.point(vh_A);
    for (auto heh : _mesh.voh_range(_mesh.from_vertex_handle(_heh)))
    {
        if (heh.to() == vh_A)
            continue;
        if (heh.next().to() == vh_A)
            continue;
        if  (heh.is_boundary())
            continue;

        auto p_B = _mesh.point(heh.to());
        auto p_C = _mesh.point(heh.next().to());

        Vec3d v_ab = p_B - p_A;
        Vec3d v_ac = p_C - p_A;
        double area = v_ab.cross(v_ac).norm() / 2;
        if  (area < _options.zero_area_delay_thres)
            return true;
    }

    return false;
}

double calc_priority(
        TriMesh& _mesh,
        HEH _heh,
        const ProgressiveMeshOptions& _options,
        OpenMesh::VPropHandleT<MatXd>* _ph_quadrics)
{
    double priority = NAN_DOUBLE;

    // Check if collapse would produce a triangle with area 0
    if (area_zero_after_collapse(_mesh, _heh, _options))
        return INF_DOUBLE;

    // Calc Priority
    switch (_options.priority)
    {
    case RANDOM:
        priority = random_priority(_heh);
        break;
    case DIHEDRAL_ANGLE_CHANGE:
        priority = dihedral_angle_change_priority(_mesh, _heh);
        break;
    case UNIFORM_DECIMATION:
        priority = uniform_priority(_mesh, _heh);
        break;
    case ERROR_QUADRICS:
        priority = error_quadrics_priority(_mesh, _heh, *_ph_quadrics);
        break;
    case ANGLE_DEFECT:
        priority = std::abs(angle_defect(_mesh, _mesh.from_vertex_handle(_heh)));
        break;
    case ANGLE_DEFECT_CHANGE_WITH_MEMORY:
        priority = angle_defect_change_with_memory_priority(_mesh, _heh);
        break;
    case ANGLE_DEFECT_CHANGE_WITHOUT_MEMORY:
        priority = angle_defect_change_without_memory_priority(_mesh, _heh);
        break;
    case MAX_MINIMAL_FACE_AREA:
        priority=  min_area_priority(_mesh, _heh);
        break;
    default:
        ISM_ERROR_throw("");
    }

    ISM_ASSERT_GEQ(priority, 0);
    ISM_ASSERT_NOT_NAN(priority);

    if (_mesh.is_boundary(_mesh.from_vertex_handle(_heh)))
        priority *= 2;

    if (_options.weight_with_minimal_face_area)
        priority += _options.weight_coeff * min_area_priority(_mesh, _heh);

    return priority;
}


struct CollapseCandidate
{
    HEH heh = HEH(-1);
    double priority = -1;

    CollapseCandidate(HEH _heh, double _priority) :
        heh(_heh),
        priority(_priority)
    { }

    bool operator<(const CollapseCandidate& _rhs) const
    {
        return priority > _rhs.priority;
    }
};

void update_queue_after_collapse_one_ring(
        TriMesh& _mesh,
        const VH _vh_to,
        const ExternalProperty<VH, bool>& _prop_locked,
        std::function<void(HEH)>& _update_priority)
{
    for (auto heh : _mesh.voh_range(_vh_to))
    {
        if (!_prop_locked[_vh_to])
            _update_priority(heh);
        if (!_prop_locked[heh.to()])
            _update_priority(heh.opp());
    }
}

void update_queue_after_collapse_two_ring(
        TriMesh& _mesh,
        const VH _vh_to,
        const ExternalProperty<VH, bool>& _prop_locked,
        std::function<void(HEH)>& _update_priority)
{
    for (auto vh : _mesh.vv_range(_vh_to))
        update_queue_after_collapse_one_ring(_mesh,
                                             vh,
                                             _prop_locked,
                                             _update_priority);
}

void update_queue_after_collapse_three_ring(
        TriMesh& _mesh,
        const VH _vh_to,
        const ExternalProperty<VH, bool>& _prop_locked,
        std::function<void(HEH)>& _update_priority)
{
    for (auto vh_1_ring :  _mesh.vv_range(_vh_to))
    {
        for (auto vh_2_ring : _mesh.vv_range(vh_1_ring))
            update_queue_after_collapse_one_ring(_mesh,
                                                 vh_2_ring,
                                                 _prop_locked,
                                                 _update_priority);
    }
}

void update_queue_after_collapse(
        TriMesh& _mesh,
        const VH _vh_to,
        const ExternalProperty<VH, bool>& _prop_locked,
        std::function<void(HEH)> _update_priority,
        const ProgressiveMeshOptions& _options)
{
    switch (_options.priority)
    {
    case RANDOM:
        break;  // Do nothing
    case DIHEDRAL_ANGLE_CHANGE:
        update_queue_after_collapse_two_ring(_mesh,
                                             _vh_to,
                                             _prop_locked,
                                             _update_priority);
        break;
    case UNIFORM_DECIMATION:
        update_queue_after_collapse_one_ring(_mesh,
                                             _vh_to,
                                             _prop_locked,
                                             _update_priority);
        break;
    case ERROR_QUADRICS:
        update_queue_after_collapse_one_ring(_mesh,
                                             _vh_to,
                                             _prop_locked,
                                             _update_priority);
        break;
    case ANGLE_DEFECT:
        update_queue_after_collapse_two_ring(_mesh,
                                             _vh_to,
                                             _prop_locked,
                                             _update_priority);
        break;
    case ANGLE_DEFECT_CHANGE_WITH_MEMORY:
        update_queue_after_collapse_three_ring(_mesh,
                                               _vh_to,
                                               _prop_locked,
                                               _update_priority);
        break;
    case ANGLE_DEFECT_CHANGE_WITHOUT_MEMORY:
        update_queue_after_collapse_three_ring(_mesh,
                                               _vh_to,
                                               _prop_locked,
                                               _update_priority);
        break;
    case MAX_MINIMAL_FACE_AREA:
        update_queue_after_collapse_one_ring(_mesh,
                                             _vh_to,
                                             _prop_locked,
                                             _update_priority);
        break;
    default:
        ISM_ERROR_throw("No valid priority");
    }
}

} // namespace

ProgressiveMesh compute_progressive_mesh_impl(
        TriMesh& _mesh,
        const int _n_target_vertices,
        const ProgressiveMeshOptions& _options,
        std::function<bool(HEH)> _is_legal_collapse,
        std::function<void(HEH)> _pre_collapse,
        std::function<void(VH, VH, VH)> _post_collapse,
        OpenMesh::VPropHandleT<MatXd>* _ph_quadrics,
        std::vector<uint>& _indep_set_sizes)
{
    TinyAD::Timer timer("Progressive Mesh");

    _mesh.request_vertex_status();
    _mesh.request_edge_status();
    _mesh.request_halfedge_status();
    _mesh.request_face_status();

    ProgressiveMesh pm;

    // Add Property to check if the priority of the queue is still valid
    OpenMesh::HPropHandleT<double> ph_current_priority;
    if (!_mesh.get_property_handle(ph_current_priority, "ph_current_priority"))
        _mesh.add_property(ph_current_priority, "ph_current_priority");

    // Initialize queue
    std::priority_queue<CollapseCandidate> queue;

    ExternalProperty<VH, bool> prop_locked(_mesh, false);

    auto update_priority = [&](const HEH _heh)
    {
        if (!_mesh.status(_heh).deleted() && _is_legal_collapse(_heh) && _mesh.is_collapse_ok(_heh))
        {
            const double property = calc_priority(_mesh, _heh, _options, _ph_quadrics);

            //#pragma omp critical
            queue.emplace(_heh, property);

            _mesh.property(ph_current_priority, _heh) = property;
        }
    };


    //#pragma omp parallel for
    for (uint i = 0; i < _mesh.n_halfedges(); ++i)
        update_priority(HEH(i));

    const int orig_num_vertices = (int)_mesh.n_vertices();
    int current_num_vertices = orig_num_vertices;

    int indep_set_size = 0;

    while (current_num_vertices > _n_target_vertices)
    {
        const CollapseCandidate cand = queue.empty() ?
                    CollapseCandidate(HEH(-1), NAN_DOUBLE)
                  : queue.top();
        if (!queue.empty())
            queue.pop();

        if (!cand.heh.is_valid()) // update queue
        {
            _indep_set_sizes.push_back(indep_set_size);
            indep_set_size = 0;

            //#pragma omp parallel for
            for (uint i = 0; i < _mesh.n_halfedges(); ++i)
                update_priority(HEH(i));

            //#pragma omp parallel for
            for(uint i = 0; i < _mesh.n_vertices(); ++i)
                prop_locked[VH(i)] = false;

            if (queue.empty())
                break;
            else
                continue;
        }

        if (prop_locked[_mesh.from_vertex_handle(cand.heh)])
            continue;

        if (_mesh.status(cand.heh).deleted())
            continue;
        if (_mesh.property(ph_current_priority, cand.heh) > cand.priority)
            continue;
        if (!_mesh.is_collapse_ok(cand.heh) || !_is_legal_collapse(cand.heh))
            continue;
        if (area_zero_after_collapse(_mesh, cand.heh, _options))
            continue;

        ISM_EXPECT(cand.priority != INF_DOUBLE);

        ISM_ASSERT(_is_legal_collapse(cand.heh));

        const auto vh_from = _mesh.from_vertex_handle(cand.heh);
        const auto vh_to = _mesh.to_vertex_handle(cand.heh);
        const auto vh_left = _mesh.opposite_vh(cand.heh); // may be invalid (boundary)
        const auto vh_right = _mesh.opposite_vh(_mesh.opposite_halfedge_handle(cand.heh)); // may be invalid (boundary)

        if (!_mesh.is_valid_handle(vh_left))
            ISM_ASSERT(_mesh.is_valid_handle(vh_right));
        if (!_mesh.is_valid_handle(vh_right))
            ISM_ASSERT(_mesh.is_valid_handle(vh_left));

        ProgressiveMesh::Record record;
        record.vh0 = vh_from;
        record.vh1 = vh_to;
        record.vhl = vh_left;
        record.vhr = vh_right;
        record.p0 = _mesh.point(vh_from);

        // Save old path ids
        _pre_collapse(cand.heh);

        if (_options.create_independent_sets)
            for (const VH vh : _mesh.vv_range(vh_from))
                prop_locked[vh] = true;

        _mesh.collapse(cand.heh);
        --current_num_vertices;

        ++indep_set_size;

        _post_collapse(vh_to, vh_left, vh_right);

        pm.log.push_back(record);

        update_queue_after_collapse(_mesh, vh_to, prop_locked, update_priority, _options);
    }

    _indep_set_sizes.push_back(indep_set_size);

    // Garbage collect the mesh while keeping track of how the vertex indices
    // are updated
    std::vector<VH> old_initial_vertices;
    for (const auto& vh : _mesh.vertices())
    {
        old_initial_vertices.push_back(vh);
    }
    std::vector<VH> new_initial_vertices = old_initial_vertices;

    std::vector<VH*> update_vh_ptr;
    std::vector<HEH*> update_heh_ptr;
    std::vector<FH*> update_fh_ptr;
    for (auto& it : new_initial_vertices)
    {
        update_vh_ptr.push_back(&it);
    }
    _mesh.garbage_collection(update_vh_ptr, update_heh_ptr, update_fh_ptr);

    // Reindexing tables
    std::vector<int> old_to_new_idx(orig_num_vertices, -1); // Size of the non-decimated mesh
    std::vector<int> new_to_old_idx(_mesh.n_vertices(), -1); // Starts at the size of the decimated mesh and then grows
    // Reindexing for vertices of the initial mesh
    for (size_t i = 0; i < old_initial_vertices.size(); ++i)
    {
        const int old_idx = old_initial_vertices[i].idx();
        const int new_idx = new_initial_vertices[i].idx();
        old_to_new_idx[old_idx] = new_idx;
        new_to_old_idx[new_idx] = old_idx;
    }
    // Reindexing for vertices of the decimation sequence
    for (auto it = pm.log.crbegin(); it != pm.log.crend(); ++it)
    {
        const auto& record = *it;
        const int old_idx = record.vh0.idx();
        const int new_idx = (int)new_to_old_idx.size();
        ISM_ASSERT(old_idx >= 0);
        ISM_ASSERT(new_idx >= 0);
        old_to_new_idx[old_idx] = new_idx;
        new_to_old_idx.push_back(old_idx);
    }
    ISM_ASSERT(old_to_new_idx.size() == new_to_old_idx.size());
    pm.orig_to_prog_idx = old_to_new_idx;
    pm.prog_to_orig_idx = new_to_old_idx;

    // Helper functions for reindexing and validation
    auto reindex = [](const VH& vh, const std::vector<int>& index_map)
    {
        if (!vh.is_valid())
            return vh;
        ISM_ASSERT(vh.idx() < (int)index_map.size());
        return VH(index_map[vh.idx()]);
    };

    auto precedes = [](const VH& vh_precursor,
                       const VH& vh_successor,
                       const std::vector<int>& index_map)
    {
        if (!vh_precursor.is_valid())
        {
            return true;
        }
        if (!vh_successor.is_valid())
        {
            return false;
        }
        return index_map[vh_precursor.idx()] <= index_map[vh_successor.idx()];
    };

    // Apply reindexing
    for (auto& record : pm.log)
    {
        ISM_ASSERT(precedes(record.vh1, record.vh0, old_to_new_idx));
        ISM_ASSERT(precedes(record.vhl, record.vh0, old_to_new_idx));
        ISM_ASSERT(precedes(record.vhr, record.vh0, old_to_new_idx));
        record.vh0 = reindex(record.vh0, old_to_new_idx);
        record.vh1 = reindex(record.vh1, old_to_new_idx);
        record.vhl = reindex(record.vhl, old_to_new_idx);
        record.vhr = reindex(record.vhr, old_to_new_idx);
    }

    return pm;
}

ProgressiveMesh compute_progressive_mesh(
        TriMesh& _mesh,
        const int _n_target_vertices,
        const ProgressiveMeshOptions& _options,
        std::vector<uint>& _indep_set_sizes,
        std::function<bool(HEH)> _is_legal_collapse,
        std::function<void(HEH)> _pre_collapse,
        std::function<void(VH, VH, VH)> _post_collapse)
{

    if (_options.priority == Priority::ERROR_QUADRICS)
    {
        ISM_DEBUG_VAR(_options.priority);
        OpenMesh::VPropHandleT<MatXd> ph_quadrics;
        if (!_mesh.get_property_handle(ph_quadrics, "PH_CURRENT_QUADRICS"))
            _mesh.add_property(ph_quadrics, "PH_CURRENT_QUADRICS");

        // Init Property
        //#pragma omp parallel for
        for (uint i = 0; i < _mesh.n_vertices(); ++i)
            _mesh.property(ph_quadrics, VH(i)) = get_quadric(_mesh, VH(i));

        auto pre_collapse_new = [&] (HEH _heh)
        {
            _mesh.property(ph_quadrics, _mesh.to_vertex_handle(_heh))
                    += _mesh.property(ph_quadrics,  _mesh.from_vertex_handle(_heh));

            _pre_collapse(_heh);
        };
        return compute_progressive_mesh_impl(_mesh, _n_target_vertices, _options, _is_legal_collapse, pre_collapse_new, _post_collapse, &ph_quadrics, _indep_set_sizes);
    }

    return compute_progressive_mesh_impl(_mesh, _n_target_vertices, _options, _is_legal_collapse, _pre_collapse, _post_collapse, nullptr, _indep_set_sizes);

}

ProgressiveMesh compute_progressive_mesh(
        TriMesh& _mesh,
        const int _n_target_vertices,
        const ProgressiveMeshOptions& _options,
        std::function<bool(HEH)> _is_legal_collapse,
        std::function<void(HEH)> _pre_collapse,
        std::function<void(VH, VH, VH)> _post_collapse)
{
    std::vector<uint> indep_set_sizes;

    return compute_progressive_mesh(_mesh, _n_target_vertices, _options, indep_set_sizes, _is_legal_collapse,
                             _pre_collapse, _post_collapse);
};



void undo_decimation_step(
        TriMesh& _mesh,
        const ProgressiveMesh::Record& _rec)
{
    ISM_ASSERT(_mesh.is_valid_handle(_rec.vh1));

    // Undo the halfedge collapse
    const auto new_heh = _mesh.vertex_split(_rec.p0, _rec.vh1, _rec.vhl, _rec.vhr);
    const auto vh0 = _mesh.from_vertex_handle(new_heh);

    ISM_ASSERT(vh0.idx() == _rec.vh0.idx());
    ISM_ASSERT(_mesh.is_manifold(_rec.vh0));
    ISM_ASSERT(_mesh.is_manifold(_rec.vh1));
}

void reconstruct_progressive_mesh(
        TriMesh& _mesh,
        const ProgressiveMesh& _pm)
{
    for (auto it = _pm.log.crbegin(); it != _pm.log.crend(); ++it)
    {
        undo_decimation_step(_mesh, *it);
    }
}

}


