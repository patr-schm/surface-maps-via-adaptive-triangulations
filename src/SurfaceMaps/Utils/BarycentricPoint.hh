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
#include <SurfaceMaps/Utils/Helpers.hh>

namespace SurfaceMaps
{

/**
 * A point on a surface expressed via barycentric coordinates in a triangle.
 * The first halfedge always points to vertex a.
 */
struct BarycentricPoint
{

BarycentricPoint() :
    heh_to_a_(-1), alpha_(NAN_DOUBLE), beta_(NAN_DOUBLE), gamma_(NAN_DOUBLE) { }

/// Applies permutation such that heh_to_a is first halfedge
BarycentricPoint(
        const HEH _heh_to_a,
        const double _alpha,
        const double _beta,
        const TriMesh& _mesh)
{
    const HEH heh_first = _mesh.halfedge_handle(_mesh.face_handle(_heh_to_a));
    if (_heh_to_a == heh_first)
    {
        alpha_ = _alpha;
        beta_ = _beta;
    }
    else if (_heh_to_a == _mesh.next_halfedge_handle(heh_first))
    {
        alpha_ = 1.0 - _alpha - _beta;
        beta_ = _alpha;
    }
    else
    {
        alpha_ = _beta;
        beta_ = 1.0 - _alpha - _beta;
    }

    heh_to_a_ = heh_first;
    gamma_ = 1.0 - alpha_ - beta_;
}

/// Applies permutation such that heh_to_a is first halfedge
BarycentricPoint(
        const HEH _heh_to_a,
        const double _alpha,
        const double _beta,
        const double _gamma,
        const TriMesh& _mesh)
{
    const HEH heh_first = _mesh.halfedge_handle(_mesh.face_handle(_heh_to_a));
    if (_heh_to_a == heh_first)
    {
        alpha_ = _alpha;
        beta_ = _beta;
        gamma_ = _gamma;
    }
    else if (_heh_to_a == _mesh.next_halfedge_handle(heh_first))
    {
        alpha_ = _gamma;
        beta_ = _alpha;
        gamma_ = _beta;
    }
    else
    {
        alpha_ = _beta;
        beta_ = _gamma;
        gamma_ = _alpha;
    }

    heh_to_a_ = heh_first;
}

BarycentricPoint(
        const FH _fh,
        const double _alpha,
        const double _beta,
        const TriMesh& _mesh) :
    heh_to_a_(_mesh.halfedge_handle(_fh)),
    alpha_(_alpha),
    beta_(_beta),
    gamma_(1.0 - _alpha - _beta)
{ }

BarycentricPoint(
        const FH _fh,
        const double _alpha,
        const double _beta,
        const double _gamma,
        const TriMesh& _mesh) :
    heh_to_a_(_mesh.halfedge_handle(_fh)),
    alpha_(_alpha),
    beta_(_beta),
    gamma_(_gamma)
{ }

BarycentricPoint(
        const Vec2d _p,
        const Vec2d _a,
        const Vec2d _b,
        const Vec2d _c,
        const HEH _heh,
        const TriMesh& _mesh)
{
    ISM_ASSERT(_heh.is_valid());
    heh_to_a_ = _heh;
    compute(_p, _a, _b, _c, heh_to_a_, _mesh, alpha_, beta_, gamma_);

    assert_first_halfedge(_mesh);
}

BarycentricPoint(
        const Vec2d _p,
        const Vec2d _a,
        const Vec2d _b,
        const Vec2d _c,
        const FH _fh,
        const TriMesh& _mesh)
{
    ISM_ASSERT(_fh.is_valid());
    heh_to_a_ = _mesh.halfedge_handle(_fh);
    compute(_p, _a, _b, _c, heh_to_a_, _mesh, alpha_, beta_, gamma_);
}

BarycentricPoint(
        const Vec2d _p,
        const FH _fh,
        const TriMesh& _mesh,
        const Parametrization& _param)
{
    ISM_ASSERT(_fh.is_valid());

    heh_to_a_ = _mesh.halfedge_handle(_fh);
    const Vec2d a = _param[vh_a(_mesh)];
    const Vec2d b = _param[vh_b(_mesh)];
    const Vec2d c = _param[vh_c(_mesh)];

    compute(_p, a, b, c, heh_to_a_, _mesh, alpha_, beta_, gamma_);
}

BarycentricPoint(
        const Vec3d _p,
        const Vec3d _a,
        const Vec3d _b,
        const Vec3d _c,
        const HEH _heh,
        const TriMesh& _mesh)
{
    ISM_ASSERT(_heh.is_valid());
    heh_to_a_ = _heh;
    compute(_p, _a, _b, _c, heh_to_a_, _mesh, alpha_, beta_, gamma_);

    assert_first_halfedge(_mesh);
}

BarycentricPoint(
        const Vec3d _p,
        const FH _fh,
        const TriMesh& _mesh)
{
    ISM_ASSERT(_fh.is_valid());

    heh_to_a_ = _mesh.halfedge_handle(_fh);
    const Vec3d a = _mesh.point(vh_a(_mesh));
    const Vec3d b = _mesh.point(vh_b(_mesh));
    const Vec3d c = _mesh.point(vh_c(_mesh));

    compute(_p, a, b, c, heh_to_a_, _mesh, alpha_, beta_, gamma_);
}

BarycentricPoint(
        const VH _vh,
        const TriMesh& _mesh)
{
    ISM_ASSERT(_vh.is_valid());

    const FH fh = *_mesh.cvf_begin(_vh);
    heh_to_a_ = _mesh.halfedge_handle(fh);

    if (_vh == _mesh.to_vertex_handle(heh_to_a_))
    {
        alpha_ = 1.0;
        beta_ = 0.0;
        gamma_ = 0.0;
    }
    else if (_vh == _mesh.to_vertex_handle(_mesh.next_halfedge_handle(heh_to_a_)))
    {
        alpha_ = 0.0;
        beta_ = 1.0;
        gamma_ = 0.0;
    }
    else if (_vh == _mesh.from_vertex_handle(heh_to_a_))
    {
        alpha_ = 0.0;
        beta_ = 0.0;
        gamma_ = 1.0;
    }
    else
        ISM_ERROR_throw("");

    assert_point(_vh, _mesh);
}

static void compute(
        const Vec2d& _p,
        const Vec2d& _a,
        const Vec2d& _b,
        const Vec2d& _c,
        const HEH _heh_to_a,
        const TriMesh& _mesh,
        double& _alpha,
        double& _beta,
        double& _gamma)
{
    const HEH heh_first = _mesh.halfedge_handle(_mesh.face_handle(_heh_to_a));
    ISM_ASSERT_EQ(_heh_to_a, heh_first);

    // Compute barycentric coordinates of 2d point
    const auto va = _a - _c;
    const auto vb = _b - _c;
    const auto vp = _p - _c;

    const double d00 = va.dot(va);
    const double d01 = va.dot(vb);
    const double d02 = va.dot(vp);
    const double d11 = vb.dot(vb);
    const double d12 = vb.dot(vp);

    const double denom = d00 * d11 - d01 * d01;

    _alpha = (d02 * d11 - d01 * d12) / denom;
    _beta = (d00 * d12 - d01 * d02) / denom;
    _gamma = 1.0 - _alpha - _beta;
}

template <typename T>
static void compute(
        const Vec2<T>& _p,
        const Vec2<T>& _a,
        const Vec2<T>& _b,
        const Vec2<T>& _c,
        const HEH _heh_to_a,
        const TriMesh& _mesh,
        double& _alpha,
        double& _beta,
        double& _gamma)
{
    const HEH heh_first = _mesh.halfedge_handle(_mesh.face_handle(_heh_to_a));
    ISM_ASSERT_EQ(_heh_to_a, heh_first);

    // Compute barycentric coordinates of 2d point
    const auto va = _a - _c;
    const auto vb = _b - _c;
    const auto vp = _p - _c;

    const double d00 = va.dot(va);
    const double d01 = va.dot(vb);
    const double d02 = va.dot(vp);
    const double d11 = vb.dot(vb);
    const double d12 = vb.dot(vp);

    const double denom = d00 * d11 - d01 * d01;

    _alpha = (d02 * d11 - d01 * d12) / denom;
    _beta = (d00 * d12 - d01 * d02) / denom;
    _gamma = 1.0 - _alpha - _beta;
}

static void compute(
        const Vec3d _p,
        const Vec3d _a,
        const Vec3d _b,
        const Vec3d _c,
        const HEH _heh_to_a,
        const TriMesh& _mesh,
        double& _alpha,
        double& _beta,
        double& _gamma)
{
    const HEH heh_first = _mesh.halfedge_handle(_mesh.face_handle(_heh_to_a));
    ISM_ASSERT_EQ(_heh_to_a, heh_first);

    // Compute barycentric coordinates of projected point
    const auto va = _a - _c;
    const auto vb = _b - _c;
    const auto vp = _p - _c;

    const double d00 = va.dot(va);
    const double d01 = va.dot(vb);
    const double d02 = va.dot(vp);
    const double d11 = vb.dot(vb);
    const double d12 = vb.dot(vp);

    const double denom = d00 * d11 - d01 * d01;

    _alpha = (d02 * d11 - d01 * d12) / denom;
    _beta = (d00 * d12 - d01 * d02) / denom;
    _gamma = 1.0 - _alpha - _beta;
}

void assert_first_halfedge(
        const TriMesh& _mesh) const
{
    const HEH heh_first = _mesh.halfedge_handle(_mesh.face_handle(heh_to_a_));
    ISM_ASSERT_EQ(heh_to_a_, heh_first);
}

Vec3d point(
        const TriMesh& _mesh) const
{
    ISM_ASSERT(is_valid());

    const Vec3d a = _mesh.point(vh_a(_mesh));
    const Vec3d b = _mesh.point(vh_b(_mesh));
    const Vec3d c = _mesh.point(vh_c(_mesh));

    return interpolate(a, b, c);
}

void assert_point(
        const VH _vh,
        const TriMesh& _mesh) const
{
    ISM_ASSERT_NORM_EPS(point(_mesh), _mesh.point(_vh), 1e-12);
}

bool is_valid() const
{
    return heh_to_a_.is_valid();
}

const double& alpha() const
{
    return alpha_;
}

const double& beta() const
{
    return beta_;
}

const double& gamma() const
{
    return gamma_;
}

bool is_inside_inclusive() const
{
    return alpha_ >= 0.0 && beta_ >= 0.0 && gamma_ >= 0.0;
}

bool is_inside_exclusive() const
{
    return alpha_ > 0.0 && beta_ > 0.0 && gamma_ > 0.0;
}

FH fh(
        const TriMesh& _mesh) const
{
    ISM_ASSERT(is_valid());
    return _mesh.face_handle(heh_to_a_);
}

void set_fh(
        const FH _fh, const TriMesh& _mesh)
{
    ISM_ASSERT(_mesh.is_valid_handle(_fh));
    set_heh(_mesh.halfedge_handle(_fh));
}

const HEH& heh() const
{
    return heh_to_a_;
}

void set_heh(
        const HEH _heh)
{
    heh_to_a_ = _heh;
}

void set_alpha_beta(
        const double& _alpha,
        const double& _beta)
{
    alpha_ = _alpha;
    beta_ = _beta;
    gamma_ = 1.0 - _alpha - _beta;
}

VH vh_a(
        const TriMesh& _mesh) const
{
    ISM_ASSERT(is_valid());
    return _mesh.to_vertex_handle(heh_to_a_);
}

VH vh_b(
        const TriMesh& _mesh) const
{
    ISM_ASSERT(is_valid());
    return _mesh.to_vertex_handle(_mesh.next_halfedge_handle(heh_to_a_));
}

VH vh_c(
        const TriMesh& _mesh) const
{
    ISM_ASSERT(is_valid());
    return _mesh.from_vertex_handle(heh_to_a_);
}

HEH heh_to_a(
        const TriMesh& _mesh) const
{
    return heh_to_a_;
}

HEH heh_to_b(
        const TriMesh& _mesh) const
{
    return _mesh.next_halfedge_handle(heh_to_a_);
}

HEH heh_to_c(
        const TriMesh& _mesh) const
{
    return _mesh.prev_halfedge_handle(heh_to_a_);
}

VH vh_closest(
        const TriMesh& _mesh) const
{
    ISM_ASSERT(is_valid());

    if (alpha_ >= beta_ && alpha_ >= gamma_)
        return vh_a(_mesh);

    if (beta_ >= alpha_ && beta_ >= gamma_)
        return vh_b(_mesh);

    return vh_c(_mesh);
}

// Returns vertex handle if exactly on vertex.
// Invalid otherwise
VH vh_on(
        const TriMesh& _mesh) const
{
    if (alpha_ == 1.0)
    {
        ISM_ASSERT_EQ(beta_, 0.0);
        ISM_ASSERT_EQ(gamma_, 0.0);
        return vh_a(_mesh);
    }
    else if (beta_ == 1.0)
    {
        ISM_ASSERT_EQ(alpha_, 0.0);
        ISM_ASSERT_EQ(gamma_, 0.0);
        return vh_b(_mesh);
    }
    else if (gamma_ == 1.0)
    {
        ISM_ASSERT_EQ(alpha_, 0.0);
        ISM_ASSERT_EQ(beta_, 0.0);
        return vh_c(_mesh);
    }
    else
        return VH(-1);
}

// Returns halfedge handle if exactly on edge.
// Invalid otherwise
// This check excludes vertices!
HEH heh_on(
        const TriMesh& _mesh) const
{
    if (alpha_ == 0.0 && beta_ != 1.0 && gamma_ != 1.0)
    {
        ISM_ASSERT_G(beta_, 0.0);
        ISM_ASSERT_L(beta_, 1.0);
        ISM_ASSERT_G(gamma_, 0.0);
        ISM_ASSERT_L(gamma_, 1.0);
        return heh_to_c(_mesh);
    }
    else if (beta_ == 0.0 && gamma_ != 1.0 && alpha_ != 1.0)
    {
        ISM_ASSERT_G(gamma_, 0.0);
        ISM_ASSERT_L(gamma_, 1.0);
        ISM_ASSERT_G(alpha_, 0.0);
        ISM_ASSERT_L(alpha_, 1.0);
        return heh_to_a(_mesh);
    }
    else if (gamma_ == 0.0 && alpha_ != 1.0 && beta_ != 1.0)
    {
        ISM_ASSERT_G(alpha_, 0.0);
        ISM_ASSERT_L(alpha_, 1.0);
        ISM_ASSERT_G(beta_, 0.0);
        ISM_ASSERT_L(beta_, 1.0);
        return heh_to_b(_mesh);
    }
    else
        return HEH(-1);
}

template <class T>
T interpolate(
        const T& _a,
        const T& _b,
        const T& _c) const
{
    if (alpha_ == 0.0)
    {
        if (beta_ == 0.0)
        {
            ISM_ASSERT_EQ(gamma_, 1.0);
            return _c;
        }
        else if (beta_ == 1.0)
        {
            ISM_ASSERT_EQ(gamma_, 0.0);
            return _b;
        }
        else
            return beta_ * _b + gamma_ * _c;
    }
    else if (beta_ == 0.0)
    {
        ISM_ASSERT_NEQ(alpha_, 0.0); // Already handled
        if (gamma_ == 0.0)
        {
            ISM_ASSERT_EQ(alpha_, 1.0);
            return _a;
        }
        else
            return alpha_ * _a + gamma_ * _c;
    }
    else if (gamma_ == 0.0)
    {
        ISM_ASSERT_NEQ(alpha_, 0.0); // Already handled
        ISM_ASSERT_NEQ(beta_, 0.0); // Already handled
        return alpha_ * _a + beta_ * _b;
    }
    else
        return alpha_ * _a + beta_ * _b + gamma_ * _c;
}

template <class T>
T interpolate(
        const ExternalProperty<VH, T>& _per_vertex,
        const TriMesh& _mesh) const
{
    ISM_ASSERT(is_valid());

    // Get  values at vertices
    const T& a = _per_vertex[vh_a(_mesh)];
    const T& b = _per_vertex[vh_b(_mesh)];
    const T& c = _per_vertex[vh_c(_mesh)];

    return interpolate(a, b, c);
}

template <class T>
T interpolate(
        const ExternalProperty<HEH, T>& _per_halfedge,
        const TriMesh& _mesh) const
{
    ISM_ASSERT(is_valid());

    // Get  values at vertices
    const T& a = _per_halfedge[heh_to_a(_mesh)];
    const T& b = _per_halfedge[heh_to_b(_mesh)];
    const T& c = _per_halfedge[heh_to_c(_mesh)];

    return interpolate(a, b, c);
}

template <class T, class VPropT>
T interpolate(
        const VPropT& _prop,
        const TriMesh& _mesh) const
{
    ISM_ASSERT(is_valid());

    // Get property values at vertices
    const T& a = _prop[vh_a(_mesh)];
    const T& b = _prop[vh_b(_mesh)];
    const T& c = _prop[vh_c(_mesh)];

    return interpolate(a, b, c);
}

private:

// This is always the first halfedge of the face.
HEH heh_to_a_;

// Store all three barycentric coordinates to faithfully represent
// points on edges and vertices.
double alpha_;
double beta_;
double gamma_;

};

}
