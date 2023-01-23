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
#include <SurfaceMaps/Viewer/IDraw.hh>
#include <glow-extras/viewer/view.hh>

namespace SurfaceMaps
{

struct GlowDraw : public IDraw
{
    GlowDraw();
    virtual ~GlowDraw() = default;

    virtual void clear() override;
    virtual void point(Vec2d _p, Color _color, const DrawStyle& _style = default_point_style) override;
    virtual void point(Vec3d _p, Color _color, const DrawStyle& _style = default_point_style) override;
    virtual void line(Vec2d _from, Vec2d _to, Color _color, const DrawStyle& _style = default_line_style) override;
    virtual void line(Vec3d _from, Vec3d _to, Color _color, const DrawStyle& _style = default_line_style) override;
    virtual void line(const HEH _heh, const TriMesh& _mesh, const Color& _color, const DrawStyle& _style = default_line_style) override;
    virtual void line(const EH _eh, const TriMesh& _mesh, const Color& _color, const DrawStyle& _style = default_line_style) override;
    virtual void triangle(Vec3d a, Vec3d b, Vec3d c, const Color& _color, const double _offset) override;
    virtual void triangle(FH _fh, const TriMesh& _mesh, const Color& _color, const double _offset) override;

    void view();

    pm::Mesh m_points;
    pm::vertex_attribute<Vec3d> pos_points;
    pm::vertex_attribute<glow::colors::color> color_points;
    pm::vertex_attribute<float> width_points;
    bool shade_points = true;
    bool world_width_points = true;
    gv::SharedPointRenderable points_renderable;
    bool points_dirty = true;

    pm::Mesh m_lines;
    pm::vertex_attribute<Vec3d> pos_lines;
    pm::edge_attribute<glow::colors::color> color_lines;
    pm::edge_attribute<float> width_lines;
    bool shade_lines = true;
    bool world_width_lines = true;
    gv::SharedLineRenderable lines_renderable;
    bool lines_dirty = true;

    pm::Mesh m_triangles;
    pm::vertex_attribute<Vec3d> pos_triangles;
    pm::face_attribute<glow::colors::color> color_triangles;
    bool shade_triangles = true;
    gv::SharedMeshRenderable triangles_renderable;
    bool triangles_dirty = true;
};

}
