/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt, Janis Born
 */

#include "MeshView.hh"
#include <SurfaceMaps/Viewer/GlowDraw.hh>
#include <SurfaceMaps/Viewer/ColorGenerator.hh>
#include <SurfaceMaps/Utils/IO.hh>
#include <SurfaceMaps/Utils/Helpers.hh>
#include <SurfaceMaps/Utils/LocalCoordinateSystem.hh>
#include <glow-extras/viewer/canvas.hh>
#include <SurfaceMaps/Viewer/HeatmapColors.hh>

namespace SurfaceMaps
{

double bb_diag(
        const TriMesh& _mesh)
{
    Vec3d min(INF_DOUBLE, INF_DOUBLE, INF_DOUBLE);
    Vec3d max(-INF_DOUBLE, -INF_DOUBLE, -INF_DOUBLE);
    for (auto v : _mesh.vertices())
    {
        min = min.cwiseMin(_mesh.point(v));
        max = max.cwiseMax(_mesh.point(v));
    }

    return (max - min).norm();
}

gv::SharedGeometricRenderable make_renderable(
        const TriMesh& _mesh)
{
    for (auto v : _mesh.vertices())
        ISM_ASSERT_FINITE_MAT(_mesh.point(v));

    // Convert to polymesh
    pm::Mesh m;
    auto pos = to_polymesh(_mesh, m);

    return gv::make_renderable(pos);
}

gv::detail::raii_config default_style()
{
    return gv::config(
                gv::no_grid,
                gv::no_outline,
                gv::no_backfacing_shadow,
                gv::maybe_empty,
                gv::background_color(tg::color3::white),
                gv::ssao_power(0.5f),
                gv::shadow_strength(0.8),
                gv::camera_fov(tg::degree(30)));
}

gv::detail::raii_config screenshot_config(
        const fs::path& _file_path,
        const glow::viewer::camera_transform& _cam_pos,
        const tg::ivec2& _size,
        const bool _transparent,
        const int _accumulation_count)
{
    make_file_directory(_file_path);

    return gv::config(
                _cam_pos,
                gv::headless_screenshot(
                               _size,
                               _accumulation_count,
                               _file_path.string(),
                               _transparent ? GL_RGBA8 : GL_RGB8));
}

gv::detail::raii_config screenshot_config(
        const fs::path& _file_path,
        const tg::ivec2& _size,
        const bool _transparent,
        const int _accumulation_count)
{
    return gv::config(
                gv::headless_screenshot(
                               _size,
                               _accumulation_count,
                               _file_path.string(),
                               _transparent ? GL_RGBA8 : GL_RGB8));
}

gv::SharedRenderable make_renderable(
        const TriMesh& _mesh,
        const glow::SharedTexture2D& _texture)
{
    ISM_ASSERT(_mesh.has_halfedge_texcoords2D());

    for (auto v : _mesh.vertices())
        ISM_ASSERT_FINITE_MAT(_mesh.point(v));

    // Convert to polymesh
    pm::Mesh m;
    auto pos = to_polymesh(_mesh, m);
    auto r = gv::make_renderable(pos);

    // Set texture coordinates
    auto uvs = m.halfedges().map([&] (auto h)
    {
        const VH vh_from(h.vertex_from().idx.value);
        const VH vh_to(h.vertex_to().idx.value);
        const HEH heh = _mesh.find_halfedge(vh_from, vh_to);
        return tg::pos2(_mesh.texcoord2D(heh)[0], _mesh.texcoord2D(heh)[1]);
    });

    // Flip texture due to different conventions
    configure(*r, gv::textured(uvs, _texture).flip());

    return r;
}

gv::SharedRenderable make_renderable(
        const TriMesh& _mesh,
        const TexCoords& _uvs,
        const glow::SharedTexture2D& _texture)
{    
    for (auto v : _mesh.vertices())
        ISM_ASSERT_FINITE_MAT(_mesh.point(v));

    // Convert to polymesh
    pm::Mesh m;
    auto pos = to_polymesh(_mesh, m);
    auto r = gv::make_renderable(pos);

    // Set texture coordinates
    auto uvs = m.halfedges().map([&] (auto h)
    {
        const VH vh_from(h.vertex_from().idx.value);
        const VH vh_to(h.vertex_to().idx.value);
        const HEH heh = _mesh.find_halfedge(vh_from, vh_to);
        return tg::pos2(_uvs[heh][0], _uvs[heh][1]);
    });

    // Flip texture due to different conventions
    configure(*r, gv::textured(uvs, _texture).flip());

    return r;
}

gv::SharedRenderable make_renderable(
        const TriMesh& _mesh,
        const ExternalProperty<FH, Color> _colors)
{
    for (auto v : _mesh.vertices())
        ISM_ASSERT_FINITE_MAT(_mesh.point(v));

    // Convert to polymesh
    pm::Mesh m;
    auto pos = to_polymesh(_mesh, m);
    auto r = gv::make_renderable(pos);

    // Set face colors
    auto tg_color = [] (auto c) { return tg::color(c[0], c[1], c[2], c[3]); };
    auto colors = m.faces().map([&] (auto f) { return tg_color(_colors[FH(f.idx.value)]); });
    configure(*r, colors);

    return r;
}

TexCoords projected_tex_coords(
        const TriMesh& _mesh,
        const Vec3d& _u_3d,
        const Vec3d& _v_3d)
{
    return TexCoords(_mesh, _mesh.halfedges().to_vector([&] (auto h)
    {
        const auto p = _mesh.point(h.to());
        return Vec2d(_u_3d.dot(p), _v_3d.dot(p));
    }));
}

void projected_tex_coords_to_property(
        TriMesh& _mesh,
        const Vec3d& _u_3d,
        const Vec3d& _v_3d)
{
    auto coords = projected_tex_coords(_mesh, _u_3d, _v_3d);

    _mesh.request_halfedge_texcoords2D();
    for (const auto h : _mesh.halfedges())
        _mesh.set_texcoord2D(h, coords[h]);
}

TexCoords normal_cluster_tex_coords(
        const TriMesh& _mesh,
        const double _scale)
{
    // Define three 2D projection planes
    std::vector planes =
    {
        std::pair<Vec3d, Vec3d> { Vec3d(0.0, 0.0, -1.0) / _scale, Vec3d(0.0, 1.0, 0.0) / _scale }, // yz plane
        std::pair<Vec3d, Vec3d> { Vec3d(1.0, 0.0, 0.0) / _scale, Vec3d(0.0, 0.0, -1.0) / _scale }, // zx plane
        std::pair<Vec3d, Vec3d> { Vec3d(1.0, 0.0, 0.0) / _scale, Vec3d(0.0, 1.0, 0.0) / _scale }, // xy plane
    };

    TexCoords uvs(_mesh);
    for (auto f : _mesh.faces())
    {
        // Classify normal (closest to x, y, or z axis)
        const Vec3d nf = _mesh.calc_face_normal(f);
        int i_max = -1;
        double dot_max = -INF_DOUBLE;
        for (int i = 0; i < 3; ++i)
        {
            const Vec3d np = planes[i].first.cross(planes[i].second).normalized();
            const double dot = std::fabs(nf.dot(np));
            if (dot > dot_max)
            {
                i_max = i;
                dot_max = dot;
            }
        }
        ISM_ASSERT_GEQ(i_max, 0);

        // Orthogonal projection
        Vec2d min_uv(INF_DOUBLE, INF_DOUBLE);
        for (auto h : f.halfedges())
        {
            const Vec3d p = _mesh.point(h.to());
            uvs[h] = Vec2d(planes[i_max].first.dot(p), planes[i_max].second.dot(p));
            min_uv = min_uv.cwiseMin(uvs[h]);
        }

        // Translate face to the r, g, or b quadrant of the texure.
        // Assumes tiling pattern:
        const double tiling_step = 1.0 / 16.0;
        auto shift = [&] (double v, double min)
        {
            // Shift to positive range
            double s = 0.0;
            if (v < min)
                s = (ceil(abs((v - min) / tiling_step))) * tiling_step;
            ISM_ASSERT_GEQ(v + s, min);

            // Left-align
            return min + std::fmod(v + s - min, tiling_step) - v;
        };

        double min_u = 0.0;
        double min_v = 0.0;
        if (i_max == 1)
        {
            min_u = 0.5;
            min_v = 0.5;
        }
        else if (i_max == 2)
        {
            min_u = 0.0;
            min_v = 0.5;
        }

        const double shift_u = shift(min_uv.x(), min_u);
        const double shift_v = shift(min_uv.y(), min_v);
        for (auto h : f.halfedges())
        {
            uvs[h] += Vec2d(shift_u, shift_v);
            ISM_ASSERT_GEQ(uvs[h].x(), 0.0);
            ISM_ASSERT_GEQ(uvs[h].y(), 0.0);
        }
    }

    return uvs;
}

void normal_cluster_tex_coords_to_property(
        TriMesh& _mesh,
        const double _scale)
{
    auto coords = normal_cluster_tex_coords(_mesh, _scale);

    _mesh.request_halfedge_texcoords2D();
    for (const auto h : _mesh.halfedges())
        _mesh.set_texcoord2D(h, coords[h]);
}

gv::detail::raii_view_closer view_mesh(
        const gv::SharedRenderable& _r,
        const std::string _caption)
{
    auto style = default_style();

    return gv::view(_r, _caption);
}

gv::detail::raii_view_closer view_mesh(
        const TriMesh& _mesh,
        const Color& _color,
        const std::string _caption)
{
    auto style = default_style();

    ExternalProperty<FH, Color> colors(_mesh, _color);
    return gv::view(make_renderable(_mesh, colors), _caption);
}

gv::SharedRenderable face_colors_renderable(
        const TriMesh& _mesh,
        const ExternalProperty<FH, Color>& _colors)
{
    // Convert to polymesh
    pm::Mesh m;
    auto pos = to_polymesh(_mesh, m);
    auto r = gv::make_renderable(pos);

    // Set face colors
    auto colors = m.faces().make_attribute<tg::color3>();
    for (auto f : m.faces())
    {
        const auto v0 = f.halfedges().first().vertex_from();
        const auto v1 = f.halfedges().first().vertex_to();
        const FH fh = _mesh.find_halfedge(VH(v0.idx.value), VH(v1.idx.value)).face();
        colors[f] = tg::color3(_colors[fh]);
    }

    return gv::make_and_configure_renderable(pos, colors);
}

void view_face_colors(
        const TriMesh& _mesh,
        const ExternalProperty<FH, Color>& _colors)
{
    auto style = default_style();
    gv::view(face_colors_renderable(_mesh, _colors));
}

void view_halfedge_colors(
        const TriMesh& _mesh,
        const ExternalProperty<HEH, Color>& _colors)
{
    // Convert to polymesh
    pm::Mesh m;
    auto pos = to_polymesh(_mesh, m);
    auto r = gv::make_renderable(pos);

    // Set corner colors
    auto colors = m.halfedges().make_attribute<tg::color3>();
    for (auto he : m.halfedges())
    {
        const auto v0 = he.vertex_from();
        const auto v1 = he.vertex_to();
        const HEH heh = _mesh.find_halfedge(VH(v0.idx.value), VH(v1.idx.value));
        colors[he] = tg::color3(_colors[heh]);
    }

    gv::view(pos, colors);
}

void draw_wireframe(
        const TriMesh& _mesh,
        const Color& _color,
        const DrawStyle& _style,
        GlowDraw& draw)
{
    for (auto e : _mesh.edges())
        draw.line(e, _mesh, _color, _style);
}

void view_wireframe(
        const TriMesh& _mesh,
        const Color& _color,
        const DrawStyle& _style)
{
    auto v = gv::view();
    auto style = default_style();

    GlowDraw draw;
    draw_wireframe(_mesh, _color, _style, draw);
    draw.view();
}

void view_embedding(
        const TriMesh& _mesh,
        const ExternalProperty<VH, Vec3d>& _pos,
        const Color& _color,
        const DrawStyle& _style)
{
    auto v = gv::view();
    auto style = default_style();

    // Convert to polymesh
    pm::Mesh m;
    auto pos = to_polymesh(_mesh, m);
    auto r = gv::make_renderable(pos);

    // Set vertex positions
    for (auto v : m.vertices())
        pos[v] = _pos[VH(v.idx.value)];

    // View mesh
    gv::view(pos);

    // View wireframe
    GlowDraw draw;
    for (auto e : _mesh.edges())
        draw.line(_pos[e.v0()], _pos[e.v1()], _color, _style);
    draw.view();
}

void view_scalar_field(
        const TriMesh& _mesh,
        const ExternalProperty<VH, double>& _field,
        const Color& _color_from,
        const Color& _color_to,
        const bool& _log_scale)
{
    auto style = default_style();

    pm::Mesh m;
    auto pos = to_polymesh(_mesh, m);

    auto v_colors = m.vertices().make_attribute<tg::color3>();
    const double min = _field.as_eigen().minCoeff();
    const double max = _field.as_eigen().maxCoeff();
    const double range = max - min;

    tg::color3 tg_color_from(_color_from[0], _color_from[1], _color_from[2]);
    tg::color3 tg_color_to(_color_to[0], _color_to[1], _color_to[2]);

    for (auto v : m.vertices())
    {
        if (range <= 0.0)
        {
            v_colors[v] = tg_color_from;
            continue;
        }

        const double val = _field[VH(v.idx.value)];

        if (_log_scale)
            v_colors[v] = tg::color3(log_color(val, min, max, _color_from, _color_to));
        else
            v_colors[v] = tg::mix(tg_color_from, tg_color_to, (val - min) / range);

    }

    gv::view(pos, v_colors);
}

void view_scalar_field(
        const TriMesh& _mesh,
        const ExternalProperty<VH, double>& _field,
        const ExternalProperty<VH, Vec3d>& _vertex_positions,
        const Color& _color_from,
        const Color& _color_to)
{
    TriMesh deformed_mesh = _mesh;
    for (auto v : _mesh.vertices())
        deformed_mesh.point(v) = _vertex_positions[v];

    view_scalar_field(deformed_mesh, _field, _color_from, _color_to);
}

void view_scalar_field(
        const TriMesh& _mesh,
        const ExternalProperty<FH, double>& _field,
        const Color& _color_from,
        const Color& _color_to)
{
    auto style = default_style();

    pm::Mesh m;
    auto pos = to_polymesh(_mesh, m);

    auto f_colors = m.faces().make_attribute<tg::color3>();
    const double min = _field.as_eigen().minCoeff();
    const double max = _field.as_eigen().maxCoeff();

    tg::color3 tg_color_from(_color_from[0], _color_from[1], _color_from[2]);
    tg::color3 tg_color_to(_color_to[0], _color_to[1], _color_to[2]);

    for (auto f : m.faces())
    {
        const double val = _field[FH(f.idx.value)];
        const double lambda = (val - min) / (max - min);
        f_colors[f] = tg::mix(tg_color_from, tg_color_to, lambda);
    }

    gv::view(pos, f_colors);
}

void view_scalar_field(
        const TriMesh& _mesh,
        const ExternalProperty<FH, double>& _field,
        const ExternalProperty<VH, Vec3d>& _vertex_positions,
        const Color& _color_from,
        const Color& _color_to)
{
    TriMesh deformed_mesh = _mesh;
    for (auto v : _mesh.vertices())
        deformed_mesh.point(v) = _vertex_positions[v];

    view_scalar_field(deformed_mesh, _field, _color_from, _color_to);
}

void view_vector_field(
        const TriMesh& _mesh,
        const ExternalProperty<VH, Vec3d>& _field,
        const ExternalProperty<VH, Vec3d>& _vertex_positions,
        const double _scale,
        const Color& _color,
        const DrawStyle& _line_style)
{
    auto style = default_style();
    auto c = gv::canvas();

    // Set color
    c.set_color(tg::color4(_color));

    // Set width
    if (_line_style.world_space)
        c.set_line_width_world(_line_style.width);
    else
        c.set_line_width_px(_line_style.width);

    for (auto v : _mesh.vertices())
        c.add_line(tg::dpos3(_vertex_positions[v]), _scale * tg::dvec3(_field[v]));
}

void view_vector_field(
        const TriMesh& _mesh,
        const ExternalProperty<VH, Vec3d>& _field,
        const double _scale,
        const Color& _color,
        const DrawStyle& _line_style)
{
    ExternalProperty<VH, Vec3d> vertex_positions(_mesh, [&] (auto v) { return _mesh.point(v); });
    view_vector_field(_mesh, _field, vertex_positions, _scale, _color, _line_style);
}

void view_vector_field(
        const TriMesh& _mesh,
        const ExternalProperty<FH, Vec3d>& _field,
        const ExternalProperty<VH, Vec3d>& _vertex_positions,
        const double _scale,
        const Color& _color,
        const DrawStyle& _line_style)
{
    auto style = default_style();
    auto c = gv::canvas();

    // Set color
    c.set_color(tg::color4(_color));

    // Set width
    if (_line_style.world_space)
        c.set_line_width_world(_line_style.width);
    else
        c.set_line_width_px(_line_style.width);

    for (auto f : _mesh.faces())
    {
        // Compute face center
        VH vh_a, vh_b, vh_c;
        handles(_mesh, f, vh_a, vh_b, vh_c);
        Vec3d p = (_vertex_positions[vh_a] + _vertex_positions[vh_b] + _vertex_positions[vh_c]) / 3.0;
        ISM_ASSERT_FINITE_MAT(p);
        ISM_ASSERT_FINITE(_scale);
        ISM_ASSERT_FINITE_MAT(p);

        c.add_line(tg::dpos3(p), tg::dvec3(_scale * _field[f]));
    }
}

void view_vector_field(
        const TriMesh& _mesh,
        const ExternalProperty<FH, Vec3d>& _field,
        const double _scale,
        const Color& _color,
        const DrawStyle& _line_style)
{
    ExternalProperty<VH, Vec3d> vertex_positions(_mesh, [&] (auto v) { return _mesh.point(v); });
    view_vector_field(_mesh, _field, vertex_positions, _scale, _color, _line_style);
}

void view_caption(const std::string& _s)
{
    gv::view(gv::make_renderable(std::vector<tg::pos3>()), gv::maybe_empty, _s);
}

void view_landmarks(
        const TriMesh& _mesh,
        const std::vector<VH>& _landmarks,
        const DrawStyle& _point_style)
{
    ColorGenerator c;

    GlowDraw draw;
    for (auto vh : _landmarks)
        draw.point(_mesh.point(vh), c.generate_next_color(), _point_style);

    draw.view();
}

void view_landmarks(
        const TriMesh& _mesh,
        const ExternalProperty<VH, Vec3d>& _embedding,
        const std::vector<VH>& _landmarks,
        const DrawStyle& _point_style)
{
    ISM_ASSERT(_embedding.size_okay(_mesh));
    ColorGenerator c;

    GlowDraw draw;
    for (auto vh : _landmarks)
        draw.point(_embedding[vh], c.generate_next_color(), _point_style);

    draw.view();
}

}
