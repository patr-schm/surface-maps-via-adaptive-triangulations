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
#include <SurfaceMaps/Utils/Filesystem.hh>
#include <SurfaceMaps/Viewer/GlowDraw.hh>
#include <SurfaceMaps/Viewer/Colors.hh>
#include <SurfaceMaps/Viewer/IDraw.hh>
#include <glow-extras/glfw/GlfwContext.hh>
#include <glow-extras/viewer/view.hh>

namespace SurfaceMaps
{

double bb_diag(
        const TriMesh& _mesh);

/// Make TriMesh renderable. Allows for gv::view(mesh).
gv::SharedGeometricRenderable make_renderable(
        const TriMesh& _mesh);

gv::detail::raii_config default_style();

gv::detail::raii_config screenshot_config(
        const fs::path& _file_path,
        const glow::viewer::camera_transform& _cam_pos,
        const tg::ivec2& _size = tg::ivec2(1920, 1080),
        const bool _transparent = true,
        const int _accumulation_count = 64);

gv::detail::raii_config screenshot_config(
        const fs::path& _file_path,
        const tg::ivec2& _size = tg::ivec2(1920, 1080),
        const bool _transparent = true,
        const int _accumulation_count = 64);

inline gv::detail::raii_config screenshot_config(
        const bool _open_viewer_instead,
        const fs::path& _file_path,
        const glow::viewer::camera_transform& _cam_pos,
        const tg::ivec2& _size = tg::ivec2(1920, 1080),
        const bool _transparent = true,
        const int _accumulation_count = 64)
{
    return _open_viewer_instead ? gv::config() : screenshot_config(_file_path, _cam_pos, _size, _transparent, _accumulation_count);
}

inline gv::detail::raii_config screenshot_config(
        const bool _open_viewer_instead,
        const fs::path& _file_path,
        const tg::ivec2& _size = tg::ivec2(1920, 1080),
        const bool _transparent = true,
        const int _accumulation_count = 64)
{
    return _open_viewer_instead ? gv::config() : screenshot_config(_file_path, _size, _transparent, _accumulation_count);
}

gv::SharedRenderable make_renderable(
        const TriMesh& _mesh,
        const glow::SharedTexture2D& _texture);

gv::SharedRenderable make_renderable(
        const TriMesh& _mesh,
        const TexCoords& _uvs,
        const glow::SharedTexture2D& _texture);

gv::SharedRenderable make_renderable(
        const TriMesh& _mesh,
        const ExternalProperty<FH, Color> _colors);

TexCoords projected_tex_coords(
        const TriMesh& _mesh,
        const Vec3d& _u_3d = Vec3d(1.0, 0.0, 0.0),
        const Vec3d& _v_3d = Vec3d(0.0, 1.0, 0.0));

void projected_tex_coords_to_property(
        TriMesh& _mesh,
        const Vec3d& _u_3d = Vec3d(1.0, 0.0, 0.0),
        const Vec3d& _v_3d = Vec3d(0.0, 1.0, 0.0));

TexCoords normal_cluster_tex_coords(
        const TriMesh& _mesh,
        const double _scale = 1.0);

void normal_cluster_tex_coords_to_property(
        TriMesh& _mesh,
        const double _scale = 1.0);

gv::detail::raii_view_closer view_mesh(
        const gv::SharedRenderable& _r,
        const std::string _caption = "");

gv::detail::raii_view_closer view_mesh(
        const TriMesh& _mesh,
        const Color& _color = WHITE,
        const std::string _caption = "");

gv::SharedRenderable face_colors_renderable(
        const TriMesh& _mesh,
        const ExternalProperty<FH, Color>& _colors);

void view_face_colors(
        const TriMesh& _mesh,
        const ExternalProperty<FH, Color>& _colors);

void view_halfedge_colors(
        const TriMesh& _mesh,
        const ExternalProperty<HEH, Color>& _colors);

void draw_wireframe(
        const TriMesh& _mesh,
        const Color& _color,
        const DrawStyle& _style,
        GlowDraw& draw);

void view_wireframe(
        const TriMesh& _mesh,
        const Color& _color = BLUE,
        const DrawStyle& = default_line_style);

void view_embedding(
        const TriMesh& _mesh,
        const ExternalProperty<VH, Vec3d>& _pos,
        const Color& _color = BLUE,
        const DrawStyle& = default_line_style);

void view_scalar_field(const TriMesh& _mesh,
        const ExternalProperty<VH, double>& _field,
        const Color& _color_from,
        const Color& _color_to,
        const bool& _log_scale = false);

void view_scalar_field(
        const TriMesh& _mesh,
        const ExternalProperty<VH, double>& _field,
        const ExternalProperty<VH, Vec3d>& _vertex_positions,
        const Color& _color_from,
        const Color& _color_to);

void view_scalar_field(
        const TriMesh& _mesh,
        const ExternalProperty<FH, double>& _field,
        const Color& _color_from,
        const Color& _color_to);

void view_scalar_field(
        const TriMesh& _mesh,
        const ExternalProperty<FH, double>& _field,
        const ExternalProperty<VH, Vec3d>& _vertex_positions,
        const Color& _color_from,
        const Color& _color_to);

void view_vector_field(
        const TriMesh& _mesh,
        const ExternalProperty<VH, Vec3d>& _field,
        const ExternalProperty<VH, Vec3d>& _vertex_positions,
        const double _scale = 1.0,
        const Color& _color = BLUE,
        const DrawStyle& _line_style = default_line_style);

void view_vector_field(
        const TriMesh& _mesh,
        const ExternalProperty<VH, Vec3d>& _field,
        const double _scale = 1.0,
        const Color& _color = BLUE,
        const DrawStyle& _line_style = default_line_style);

void view_vector_field(
        const TriMesh& _mesh,
        const ExternalProperty<FH, Vec3d>& _field,
        const ExternalProperty<VH, Vec3d>& _vertex_positions,
        const double _scale = 1.0,
        const Color& _color = BLUE,
        const DrawStyle& _line_style = default_line_style);

void view_vector_field(
        const TriMesh& _mesh,
        const ExternalProperty<FH, Vec3d>& _field,
        const double _scale = 1.0,
        const Color& _color = BLUE,
        const DrawStyle& _line_style = default_line_style);

void view_caption(
        const std::string& _s);

void view_landmarks(
        const TriMesh& _mesh,
        const std::vector<VH>& _landmarks,
        const DrawStyle& _point_style = default_point_style);

void view_landmarks(
        const TriMesh& _mesh,
        const ExternalProperty<VH, Vec3d>& _embedding,
        const std::vector<VH>& _landmarks,
        const DrawStyle& _point_style = default_point_style);

}
