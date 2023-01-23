/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt, Janis Born
 */

#include "GlowDraw.hh"

#include <SurfaceMaps/Utils/Helpers.hh>

namespace SurfaceMaps
{

GlowDraw::GlowDraw()
{
    pos_points = m_points.vertices().make_attribute<Vec3d>();
    color_points = m_points.vertices().make_attribute<glow::colors::color>();
    width_points = m_points.vertices().make_attribute<float>();

    pos_lines = m_lines.vertices().make_attribute<Vec3d>();
    color_lines = m_lines.edges().make_attribute<glow::colors::color>();
    width_lines = m_lines.edges().make_attribute<float>();

    pos_triangles = m_triangles.vertices().make_attribute<Vec3d>();
    color_triangles = m_triangles.faces().make_attribute<glow::colors::color>();
}

void GlowDraw::clear()
{
    m_points.clear();
    pos_points.clear();
    color_points.clear();
    width_points.clear();
    shade_points = true;
    world_width_points = true;
    points_renderable.reset();
    points_dirty = true;

    m_lines.clear();
    pos_lines.clear();
    color_lines.clear();
    width_lines.clear();
    shade_lines = true;
    world_width_lines = true;
    lines_renderable.reset();
    lines_dirty = true;

    m_triangles.clear();
    pos_triangles.clear();
    color_triangles.clear();
    shade_triangles = true;
    triangles_renderable.reset();
    triangles_dirty = true;
}

void GlowDraw::point(Vec2d _p, Color _color, const DrawStyle& _style)
{
    point(Vec3d(_p[0], _p[1], 0.0), _color, _style);
}

void GlowDraw::point(Vec3d _p, Color _color, const DrawStyle& _style)
{
    auto v = m_points.vertices().add();
    pos_points[v] = _p;
    color_points[v] = glow::colors::color(_color[0], _color[1], _color[2], _color[3]);
    width_points[v] = _style.width;
    world_width_points = _style.world_space;
    shade_points = _style.shaded;
    points_dirty = true;
}

void GlowDraw::line(Vec2d _from, Vec2d _to, Color _color, const DrawStyle& _style)
{
    line(Vec3d(_from[0], _from[1], 0.0), Vec3d(_to[0], _to[1], 0.0), _color, _style);
}

void GlowDraw::line(Vec3d _from, Vec3d _to, Color _color, const DrawStyle& _style)
{
    auto v_from = m_lines.vertices().add();
    auto v_to = m_lines.vertices().add();
    pos_lines[v_from] = _from;
    pos_lines[v_to] = _to;
    auto e = m_lines.edges().add_or_get(v_from, v_to);
    color_lines[e] = glow::colors::color(_color[0], _color[1], _color[2], _color[3]);
    width_lines[e] = _style.width;
    world_width_lines = _style.world_space;
    shade_lines = _style.shaded;
    lines_dirty = true;
}

void GlowDraw::line(const HEH _heh, const TriMesh& _mesh, const Color& _color, const DrawStyle& _style)
{
    auto h = make_smart(_heh, _mesh);
    line(_mesh.point(h.from()), _mesh.point(h.to()), _color, _style);
}

void GlowDraw::line(const EH _eh, const TriMesh& _mesh, const Color& _color, const DrawStyle& _style)
{
    line(_mesh.halfedge_handle(_eh, 0), _mesh, _color, _style);
}

void GlowDraw::triangle(Vec3d a, Vec3d b, Vec3d c, const Color& _color, const double _offset)
{
    // Offset points in normal direction
    const auto n = (b-a).cross(c-a).normalized();
    a += _offset * n;
    b += _offset * n;
    c += _offset * n;

    const auto va = m_triangles.vertices().add();
    const auto vb = m_triangles.vertices().add();
    const auto vc = m_triangles.vertices().add();
    const auto f = m_triangles.faces().add(va, vb, vc);
    pos_triangles[va] = a;
    pos_triangles[vb] = b;
    pos_triangles[vc] = c;
    color_triangles[f] = glow::colors::color(_color[0], _color[1], _color[2], _color[3]);
    shade_triangles = true;
    triangles_dirty = true;
}

void GlowDraw::triangle(FH _fh, const TriMesh& _mesh, const Color& _color, const double _offset)
{
    VH va, vb, vc;
    handles(_mesh, _fh, va, vb, vc);
    triangle(_mesh.point(va), _mesh.point(vb), _mesh.point(vc), _color, _offset);
}

void GlowDraw::view()
{
    // Points
    if (points_dirty)
    {
        if (world_width_points)
            points_renderable = gv::make_renderable(gv::points(pos_points).point_size_world(width_points));
        else
            points_renderable = gv::make_renderable(gv::points(pos_points).point_size_px(width_points));
        configure(*points_renderable, color_points);
        points_dirty = false;
    }

    if (points_renderable)
    {
        if (shade_points)
            gv::view(points_renderable, gv::maybe_empty);
        else
            gv::view(points_renderable, gv::maybe_empty, gv::no_shading);
    }

    // Lines
    if (lines_dirty)
    {
        if (world_width_lines)
            lines_renderable = gv::make_renderable(gv::lines(pos_lines).line_width_world(width_lines));
        else
            lines_renderable = gv::make_renderable(gv::lines(pos_lines).line_width_px(width_lines));
        configure(*lines_renderable, color_lines);
        lines_dirty = false;
    }

    if (lines_renderable)
    {
        if (shade_lines)
            gv::view(lines_renderable, gv::maybe_empty);
        else
            gv::view(lines_renderable, gv::maybe_empty, gv::no_shading);
    }

    // Triangles
    if (triangles_dirty)
    {
        triangles_renderable = gv::make_renderable(gv::polygons(pos_triangles));
        configure(*triangles_renderable, color_triangles);
        triangles_dirty = false;
    }

    if (triangles_renderable)
    {
        if (shade_triangles)
            gv::view(triangles_renderable, gv::maybe_empty);
        else
            gv::view(triangles_renderable, gv::maybe_empty, gv::no_shading);
    }
}

}
