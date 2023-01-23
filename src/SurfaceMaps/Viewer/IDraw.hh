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
#include <SurfaceMaps/Viewer/Colors.hh>

namespace SurfaceMaps
{

struct DrawStyle
{
    bool world_space = true; // vs screen space
    float width = 0.0025;
    bool shaded = true;
};

struct WidthWorld : public DrawStyle
{
    explicit WidthWorld(
            const float _width)
    {
        world_space = true;
        width = _width;
    }

    WidthWorld(
            const float _width,
            const bool _shaded)
    {
        world_space = true;
        width = _width;
        shaded = _shaded;
    }
};

struct WidthScreen : public DrawStyle
{
    explicit WidthScreen(
            const float _width)
    {
        world_space = false;
        width = _width;
    }

    WidthScreen(
            const float _width,
            const bool _shaded)
    {
        world_space = false;
        width = _width;
        shaded = _shaded;
    }
};

static const WidthScreen default_line_style(0.5, true);
static const WidthScreen default_point_style(4, true);

class IDraw
{

public:
    IDraw() = default;
    virtual ~IDraw() = default;

    virtual void clear() = 0;
    virtual void point(Vec2d _p, Color _color, const DrawStyle& _style = default_point_style) = 0;
    virtual void point(Vec3d _p, Color _color, const DrawStyle& _style = default_point_style) = 0;
    virtual void line(Vec2d _from, Vec2d _to, Color _color, const DrawStyle& _style = default_line_style) = 0;
    virtual void line(Vec3d _from, Vec3d _to, Color _color, const DrawStyle& _style = default_line_style) = 0;
    virtual void line(const HEH _heh, const TriMesh& _mesh, const Color& _color, const DrawStyle& _style = default_line_style) = 0;
    virtual void line(const EH _eh, const TriMesh& _mesh, const Color& _color, const DrawStyle& _style = default_line_style) = 0;
    virtual void triangle(Vec3d a, Vec3d b, Vec3d c, const Color& _color, const double _offset) = 0;
    virtual void triangle(FH _fh, const TriMesh& _mesh, const Color& _color, const double _offset) = 0;
};

class DummyDraw : public IDraw
{
public:
    DummyDraw() = default;
    virtual ~DummyDraw() = default;

    virtual void clear() override { }
    virtual void point(Vec2d _p, Color _color, const DrawStyle& _style = default_point_style) override { }
    virtual void point(Vec3d _p, Color _color, const DrawStyle& _style = default_point_style) override { }
    virtual void line(Vec2d _from, Vec2d _to, Color _color, const DrawStyle& _style = default_line_style) override { }
    virtual void line(Vec3d _from, Vec3d _to, Color _color, const DrawStyle& _style = default_line_style) override { }
    virtual void line(const HEH _heh, const TriMesh& _mesh, const Color& _color, const DrawStyle& _style = default_line_style) override { }
    virtual void line(const EH _eh, const TriMesh& _mesh, const Color& _color, const DrawStyle& _style = default_line_style) override { }
    virtual void triangle(Vec3d a, Vec3d b, Vec3d c, const Color& _color, const double _offset) override { }
    virtual void triangle(FH _fh, const TriMesh& _mesh, const Color& _color, const double _offset) override { }
};

}
