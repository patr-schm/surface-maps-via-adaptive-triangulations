/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Dörte Pieper, Patrick Schmidt
 */

#include <SurfaceMaps/AdaptiveTriangulations/Visualization.hh>

#include <SurfaceMaps/Viewer/MeshView.hh>
#include <SurfaceMaps/Viewer/ConstantCurvature.hh>
#include <SurfaceMaps/AdaptiveTriangulations/Helpers.hh>
#include <SurfaceMaps/AdaptiveTriangulations/DistortionHeatmap.hh>
#include <SurfaceMaps/AdaptiveTriangulations/LiftToSurface.hh>
#include <SurfaceMaps/Viewer/ColorGenerator.hh>
#include <SurfaceMaps/Utils/IO.hh>
#include <glow-extras/viewer/canvas.hh>

namespace SurfaceMaps
{

void view_texture_frontal_projection_T(
        const TriMesh& _mesh_T_A,
        const TriMesh& _mesh_T_B,
        const std::string& _name,
        const int &_projection_dir,
        const double& _texture_factor,
        const fs::path& _texture_path,
        const std::experimental::filesystem::__cxx11::path &screenshot_path,
        const tg::ivec2 &_size,
        const bool _transparent,
        const glow::viewer::camera_transform &_cam_pos)
{
    ExternalProperty<HEH, Vec2d> texcoords_T(_mesh_T_A);
    for (auto heh : _mesh_T_A.halfedges()) // leave z coord
    {
        texcoords_T[heh] = Vec2d(_mesh_T_A.point(heh.to())[(_projection_dir+1) % 3]*_texture_factor, _mesh_T_A.point(heh.to())[(_projection_dir+2) % 3]*_texture_factor);
    }
    if (!screenshot_path.empty())
        fs::create_directories(screenshot_path.parent_path());
    auto s = screenshot_path.empty() ? gv::config() : screenshot_config(screenshot_path, _cam_pos,_size,_transparent);

    glow::SharedTexture2D texture = read_texture(_texture_path);
    auto g = gv::grid();
    {
        auto v = gv::view();
        view_caption(_name);
        view_mesh(make_renderable(_mesh_T_A, texcoords_T, texture));
    }
    {
        auto v = gv::view();
        view_caption(_name);
        view_mesh(make_renderable(_mesh_T_B, texcoords_T, texture));
    }

}

void view_texture_frontal_projection_T(
        const MapState& _map_state,
        const int _ref_mesh_num,
        const std::vector<int>& _indices,
        const int &_projection_dir,
        const double& _texture_factor,
        const fs::path& _texture_path,
        const std::experimental::filesystem::__cxx11::path & _screenshot_path,
        const tg::ivec2 &_size,
        const bool _transparent,
        const glow::viewer::camera_transform &_cam_pos)
{
    if (_ref_mesh_num >= (int) _map_state.meshes_input.size())
    {
        ISM_WARNING("Desired number for frontal projection too large. Use 0 instead");
        view_texture_frontal_projection_T(_map_state, 0, _indices, _projection_dir, _texture_factor, _texture_path, _screenshot_path, _size, _transparent, _cam_pos);
        return;
    }

    // Compute texture coordinates for first T mesh to map it onto others
    TriMesh lifted_T_ref = embedding_to_mesh(_map_state.mesh_T, lift_to_surface(_map_state.meshes_input[_ref_mesh_num], _map_state.mesh_T, _map_state.embeddings_input[_ref_mesh_num], _map_state.embeddings_T[_ref_mesh_num], _map_state.bsp_embeddings_input[_ref_mesh_num]));
    ExternalProperty<HEH, Vec2d> texcoords_T(_map_state.mesh_T);
    for (auto heh : _map_state.mesh_T.halfedges()) // leave z coord
    {
        texcoords_T[heh] = Vec2d(lifted_T_ref.point(heh.to())[(_projection_dir+1) % 3]*_texture_factor, lifted_T_ref.point(heh.to())[(_projection_dir+2) % 3]*_texture_factor);
    }

    glow::SharedTexture2D texture = read_texture(_texture_path);
    {
        for(int i= 0; i< (int)_map_state.meshes_input.size(); i++)
        {
            if (!_screenshot_path.empty())
                fs::create_directories(_screenshot_path.parent_path());

            if(_indices.size() == _map_state.meshes_input.size())
            {
                auto s = _screenshot_path.empty() ? gv::config() : screenshot_config(_screenshot_path.parent_path() / (_screenshot_path.stem().string() + ("_" + std::to_string(_indices[i]) + ".png")), _cam_pos,_size,_transparent);
                TriMesh lifted_T_i = embedding_to_mesh(_map_state.mesh_T, lift_to_surface(_map_state.meshes_input[i], _map_state.mesh_T, _map_state.embeddings_input[i], _map_state.embeddings_T[i], _map_state.bsp_embeddings_input[i]));
                auto v = gv::view();
                view_mesh(make_renderable(lifted_T_i, texcoords_T, texture));
            }
            else
            {
                auto s = _screenshot_path.empty() ? gv::config() : screenshot_config(_screenshot_path.parent_path() / (_screenshot_path.stem().string() + ("_" + std::to_string(i) + ".png")), _cam_pos,_size,_transparent);
                TriMesh lifted_T_i = embedding_to_mesh(_map_state.mesh_T, lift_to_surface(_map_state.meshes_input[i], _map_state.mesh_T, _map_state.embeddings_input[i], _map_state.embeddings_T[i], _map_state.bsp_embeddings_input[i]));
                auto v = gv::view();
                view_mesh(make_renderable(lifted_T_i, texcoords_T, texture));
            }
        }
    }
}

/// Show texture mapped to target mesh
void view_texture_frontal_projection_input(
        const MapState& _map_state,
        const int _source_mesh_idx,
        const int _target_mesh_idx,
        const int &_projection_dir,
        const double& _texture_factor,
        const glow::SharedTexture2D& _texture)
{
    // Compute texture coordinates on source mesh
    ExternalProperty<HEH, Vec2d> texcoords_source(_map_state.meshes_input[_source_mesh_idx]);
    for (auto heh : _map_state.meshes_input[_source_mesh_idx].halfedges())
    {
        texcoords_source[heh] = Vec2d(
                    _map_state.meshes_input[_source_mesh_idx].point(heh.to())[(_projection_dir+1) % 3]*_texture_factor,
                    _map_state.meshes_input[_source_mesh_idx].point(heh.to())[(_projection_dir+2) % 3]*_texture_factor);
    }

    // Transfer texture coordinates to target mesh
    ExternalProperty<HEH, Vec2d> texcoords_target(_map_state.meshes_input[_target_mesh_idx]);
    TriMesh mesh_embedding_T_target = embedding_to_mesh(_map_state.mesh_T, _map_state.embeddings_T[_target_mesh_idx]);
    BSPTree bsp_T_i(mesh_embedding_T_target);
    for (auto heh : _map_state.meshes_input[_target_mesh_idx].halfedges())
    {
        // Compute barycentric coordinates of points of target mesh in T
        SVH v_to = heh.to();
        SFH fh_T;
        double alpha_T, beta_T, gamma_T;
        bsp_tree_barys_face(_map_state.embeddings_input[_target_mesh_idx][v_to], mesh_embedding_T_target, bsp_T_i, alpha_T, beta_T, gamma_T, fh_T);
        BarycentricPoint bary_T(fh_T, alpha_T, beta_T, _map_state.mesh_T);

        // Map points of mesh i to ref_mesh sphere
        const Vec3d p_emb_on_ref = bary_T.interpolate(_map_state.embeddings_T[_source_mesh_idx], _map_state.mesh_T).normalized();

        // Compute barycentric coordinates of target mesh in ref_mesh triangle
        SFH fh_ref;
        double alpha_ref, beta_ref, gamma_ref;
        bsp_tree_barys_face(p_emb_on_ref, _map_state.meshes_embeddings_input[_source_mesh_idx], _map_state.bsp_embeddings_input[_source_mesh_idx], alpha_ref, beta_ref, gamma_ref, fh_ref);

        // Compute corresponding texcoords
        SHEH heh_a, heh_b, heh_c;
        handles(_map_state.meshes_input[_source_mesh_idx], fh_ref, heh_a, heh_b, heh_c);
        texcoords_target[heh] = texcoords_source[heh_a] * alpha_ref + texcoords_source[heh_b] * beta_ref + texcoords_source[heh_c] * gamma_ref;
    }

    view_mesh(make_renderable(_map_state.meshes_input[_target_mesh_idx], texcoords_target, _texture));
}

void view_texture_frontal_projection_input(
        const MapState& _map_state,
        const int _ref_mesh_num,
        const std::vector<int>& _indices,
        const int &_projection_dir,
        const double& _texture_factor,
        const fs::path& _texture_path,
        const std::experimental::filesystem::__cxx11::path & _screenshot_path,
        const tg::ivec2 &_size,
        const bool _transparent,
        const glow::viewer::camera_transform &_cam_pos)
{
    if (_ref_mesh_num >= (int) _map_state.meshes_input.size())
    {
        ISM_WARNING("Desired number for frontal projection too large. Use 0 instead");
        view_texture_frontal_projection_input(_map_state, 0, _indices, _projection_dir, _texture_factor, _texture_path, _screenshot_path, _size, _transparent, _cam_pos);
        return;
    }

    // Compute texture coordinates for first T mesh to map it onto others
    ExternalProperty<HEH, Vec2d> texcoords_ref(_map_state.meshes_input[_ref_mesh_num]);
    for (auto heh : _map_state.meshes_input[_ref_mesh_num].halfedges()) // leave z coord
    {
        texcoords_ref[heh] = Vec2d(_map_state.meshes_input[_ref_mesh_num].point(heh.to())[(_projection_dir+1) % 3]*_texture_factor, _map_state.meshes_input[_ref_mesh_num].point(heh.to())[(_projection_dir+2) % 3]*_texture_factor);
    }

    glow::SharedTexture2D texture = read_texture(_texture_path);
    {
        for(int i= 0; i< (int)_map_state.meshes_input.size(); i++)
        {
            ExternalProperty<HEH, Vec2d> texcoords_i(_map_state.meshes_input[i]);
            if (i == _ref_mesh_num)
            {
                texcoords_i = texcoords_ref;
            }
            else
            {
                TriMesh mesh_embedding_T_i = embedding_to_mesh(_map_state.mesh_T, _map_state.embeddings_T[i]);
                BSPTree bsp_T_i(mesh_embedding_T_i);

                for (auto heh : _map_state.meshes_input[i].halfedges())
                {
                    // Compute barycentric coordinates of points of mesh i in T
                    SVH v_to = heh.to();
                    SFH fh_T;
                    double alpha_T, beta_T, gamma_T;
                    bsp_tree_barys_face(_map_state.embeddings_input[i][v_to], mesh_embedding_T_i, bsp_T_i, alpha_T, beta_T, gamma_T, fh_T);
                    BarycentricPoint bary_T(fh_T, alpha_T, beta_T, _map_state.mesh_T);

                    // Map points of mesh i to ref_mesh sphere
                    const Vec3d p_emb_on_ref = bary_T.interpolate(_map_state.embeddings_T[_ref_mesh_num], _map_state.mesh_T).normalized();

                    // Compute barycentric coordinates of mesh i in ref_mesh triangle
                    SFH fh_ref;
                    double alpha_ref, beta_ref, gamma_ref;
                    bsp_tree_barys_face(p_emb_on_ref, _map_state.meshes_embeddings_input[_ref_mesh_num], _map_state.bsp_embeddings_input[_ref_mesh_num], alpha_ref, beta_ref, gamma_ref, fh_ref);

                    // Compute corresponding texcoords
                    SHEH heh_a, heh_b, heh_c;
                    handles(_map_state.meshes_input[_ref_mesh_num], fh_ref, heh_a, heh_b, heh_c);
                    texcoords_i[heh] = texcoords_ref[heh_a] * alpha_ref + texcoords_ref[heh_b] * beta_ref + texcoords_ref[heh_c] * gamma_ref;
                }
            }

            if (!_screenshot_path.empty())
                fs::create_directories(_screenshot_path.parent_path());

            if(_indices.size() == _map_state.meshes_input.size())
            {
                auto s = _screenshot_path.empty() ? gv::config() : screenshot_config(_screenshot_path.parent_path() / (_screenshot_path.stem().string() + ("_" + std::to_string(_indices[i]) + ".png")), _cam_pos,_size,_transparent);
                auto v = gv::view();
                view_mesh(make_renderable(_map_state.meshes_input[i], texcoords_i, texture));
            }
            else
            {
                auto s = _screenshot_path.empty() ? gv::config() : screenshot_config(_screenshot_path.parent_path() / (_screenshot_path.stem().string() + ("_" + std::to_string(i) + ".png")), _cam_pos,_size,_transparent);
                auto v = gv::view();
                view_mesh(make_renderable(_map_state.meshes_input[i], texcoords_i, texture));
            }
        }
    }
}

/// Show texture mapped to target mesh
void view_texture_halfedgetexcoords_input(
        const MapState& _map_state,
        const int _source_mesh_idx,
        const int _target_mesh_idx,
        const glow::SharedTexture2D& _texture)
{
    ISM_ASSERT(_map_state.meshes_input[_source_mesh_idx].has_halfedge_texcoords2D());

    // Transfer texture coordinates to target mesh
    ExternalProperty<HEH, Vec2d> texcoords_target(_map_state.meshes_input[_target_mesh_idx]);
    TriMesh mesh_embedding_T_target = embedding_to_mesh(_map_state.mesh_T, _map_state.embeddings_T[_target_mesh_idx]);
    BSPTree bsp_T_i(mesh_embedding_T_target);
    for (auto heh : _map_state.meshes_input[_target_mesh_idx].halfedges())
    {
        // Compute barycentric coordinates of points of target mesh in T
        SVH v_to = heh.to();
        SFH fh_T;
        double alpha_T, beta_T, gamma_T;
        bsp_tree_barys_face(_map_state.embeddings_input[_target_mesh_idx][v_to], mesh_embedding_T_target, bsp_T_i, alpha_T, beta_T, gamma_T, fh_T);
        BarycentricPoint bary_T(fh_T, alpha_T, beta_T, _map_state.mesh_T);

        // Map points of mesh i to ref_mesh sphere
        const Vec3d p_emb_on_ref = bary_T.interpolate(_map_state.embeddings_T[_source_mesh_idx], _map_state.mesh_T).normalized();

        // Compute barycentric coordinates of target mesh in ref_mesh triangle
        SFH fh_ref;
        double alpha_ref, beta_ref, gamma_ref;
        bsp_tree_barys_face(p_emb_on_ref, _map_state.meshes_embeddings_input[_source_mesh_idx], _map_state.bsp_embeddings_input[_source_mesh_idx], alpha_ref, beta_ref, gamma_ref, fh_ref);

        // Compute corresponding texcoords
        SHEH heh_a, heh_b, heh_c;
        handles(_map_state.meshes_input[_source_mesh_idx], fh_ref, heh_a, heh_b, heh_c);
        texcoords_target[heh] = _map_state.meshes_input[_source_mesh_idx].texcoord2D(heh_a) * alpha_ref
                + _map_state.meshes_input[_source_mesh_idx].texcoord2D(heh_b) * beta_ref
                + _map_state.meshes_input[_source_mesh_idx].texcoord2D(heh_c) * gamma_ref;
    }

    view_mesh(make_renderable(_map_state.meshes_input[_target_mesh_idx], texcoords_target, _texture));
}

void view_lifted_with_landmarks_in_meshes(
        const MapState& _map_state,
        const fs::path& _screenshot_path,
        const tg::ivec2& _size,
        const bool _transparent,
        const glow::viewer::camera_transform& _cam_pos)
{
    if (!_screenshot_path.empty())
        fs::create_directories(_screenshot_path.parent_path());
    auto s = _screenshot_path.empty() ? gv::config() : screenshot_config(_screenshot_path, _cam_pos,_size,_transparent);

    float world_size = 0.005;
    gv::arrow_style arrow_style;
    // Minimal arrow length with default style, if arrow length is shorter not correct visualiation because "not enough space: start from to and go backwards" in canvas.cc
    auto min_arrow_length = world_size * (arrow_style.length_factor + arrow_style.shaft_min_length_factor + arrow_style.margin_arrow_factor + arrow_style.margin_shaft_factor);

    auto g = gv::grid();
    {
        for(int i= 0; i< (int)_map_state.meshes_input.size(); i++)
        {
            TriMesh lifted_T_i = embedding_to_mesh(_map_state.mesh_T, lift_to_surface(_map_state.meshes_input[i], _map_state.mesh_T, _map_state.embeddings_input[i], _map_state.embeddings_T[i], _map_state.bsp_embeddings_input[i]));
            {
                auto v = gv::view();
                auto c = gv::canvas();
                view_mesh(lifted_T_i, Color(1.0, 1.0, 1.0, 1.0));
                view_mesh(_map_state.meshes_input[i], Color(1.0, 1.0, 1.0, 0.3));
                view_wireframe(lifted_T_i, MAGENTA);
                //view_wireframe(_map_state.meshes_input[i], Color(1.0, 1.0, 1.0, 0.7));
                ColorGenerator colors;
                for (auto j=0; j< (int)_map_state.landmarks_input[i].size(); j++)
                {
                    auto vh_T = _map_state.landmarks_T[j];
                    auto vh_A = _map_state.landmarks_input[i][j];
                    auto col = colors.generate_next_color();
                    c.add_point(lifted_T_i.point(vh_T), tg::color4(col));
                    auto test = tg::pos3(lifted_T_i.point(vh_T)) - tg::pos3(_map_state.meshes_input[i].point(vh_A));
                    if (length(test) >= min_arrow_length)
                        c.add_arrow(tg::pos3(lifted_T_i.point(vh_T)), tg::pos3(_map_state.meshes_input[i].point(vh_A)), world_size, tg::color3(col));
                    // Arrow has to be drawn thinner etc. as it would otherwise not fit for the length
                    else if (length(test) > 0.0)
                    {
                        c.add_arrow(tg::pos3(lifted_T_i.point(vh_T)), tg::pos3(_map_state.meshes_input[i].point(vh_A)), world_size * (length(test)/min_arrow_length), tg::color3(col));
                    }
                }
            }
        }
    }
}


void view_landmarks_on_spheres(
        const MapState& _map_state,
        const fs::path& _screenshot_path,
        const tg::ivec2& _size,
        const bool _transparent,
        const glow::viewer::camera_transform& _cam_pos)
{
    if (!_screenshot_path.empty())
        fs::create_directories(_screenshot_path.parent_path());
    auto s = _screenshot_path.empty() ? gv::config() : screenshot_config(_screenshot_path, _cam_pos,_size,_transparent);

    DrawStyle point_style = default_point_style;

    auto g = gv::grid();

    for (int i = 0; i < (int)_map_state.meshes_input.size(); ++i)
    {
        auto v = gv::view();
        view_ccm_embedding(_map_state.meshes_input[i], _map_state.embeddings_input[i], Spherical, BLUE_50);
        view_ccm_embedding(_map_state.mesh_T, _map_state.embeddings_T[i], Spherical, MAGENTA);
        point_style.width = 8.0;
        view_landmarks(_map_state.meshes_input[i], _map_state.embeddings_input[i], _map_state.landmarks_input[i], point_style);
        point_style.width = 16.0;
        view_landmarks(_map_state.mesh_T, _map_state.embeddings_T[i], _map_state.landmarks_T, point_style);
    }
}

}
