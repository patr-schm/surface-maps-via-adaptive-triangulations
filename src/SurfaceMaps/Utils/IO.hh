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

namespace polymesh
{

class Mesh;
struct vertex_handle;
struct edge_handle;
template <class AttrT> struct vertex_attribute;

}
namespace pm = polymesh;

namespace tg
{

template <int d, typename T> struct pos;
using pos3 = pos<3, float>;

}

namespace glow
{
    class Texture2D;
}

namespace SurfaceMaps
{

/// Create directory (recursive)
void make_directory(
        const std::string& _dir_path);

void make_file_directory(
        const std::string& _file_path);

std::ifstream open_file_read(
        const std::string& _path,
        const bool _throw_if_not_exist = true);

void remove_comment(
        std::string& str,
        std::string comment);

void trim(
        std::string& str);

std::vector<std::string> read_line(
        std::ifstream& _file,
        const char _delimiter);

void append_line(
        const fs::path& _file_path,
        const std::string& _line,
        const bool _silent = false);

template<typename Arg0, typename ...Args>
void append_to_csv(
        const fs::path& _file_path,
        Arg0 arg0,
        Args ...args)
{
    std::stringstream s;
    s << arg0;
    ((s << ", " << args), ...);
    append_line(_file_path, s.str());
}

template<typename Arg0, typename ...Args>
void append_to_csv_silent(
        const fs::path& _file_path,
        Arg0 arg0,
        Args ...args)
{
    std::stringstream s;
    s << arg0;
    ((s << ", " << args), ...);
    append_line(_file_path, s.str(), true);
}

/// Read OpenMesh with per-corner texture coordinates
TriMesh read_mesh(
        const fs::path& _path,
        const bool _silent = false);

/// Write OpenMesh with per-corner texture coordinates
void write_mesh(
        const TriMesh& _mesh,
        const fs::path& _file_path);

/// Get per-corner texture coordinates from mesh
TexCoords tex_coords(
        const TriMesh& _mesh);

enum class LandmarkFormat
{
    Pinned,  // Default: A list of vertex indices
    VTS,     // The format used by the annotated SHREC07 dataset
    Unknown, // Unrecognized format
};

enum class LandmarkType
{
    Init,   // li: Use only for map init, discard afterwards
    PreOpt, // lp: Use for map init and pre-optimization, discard afterwards
    Keep,   // lk: Use in init, pre-optimization and optimization
};

LandmarkFormat guess_landmark_format(const fs::path& _path);

/// Read landmarks from file.
/// Tries to guess the landmark format from the file extension.
std::vector<VH> read_landmarks(
        const fs::path& _path,
        const LandmarkType _type = LandmarkType::Keep,
        const bool _throw_if_not_exist = true);

std::vector<VH> read_landmarks_pinned(
        const fs::path& _path,
        const LandmarkType _type = LandmarkType::Keep,
        const bool _throw_if_not_exist = true);
std::vector<VH> read_landmarks_vts(
        const fs::path& _path,
        const bool _throw_if_not_exist = true);

/// Write landmarks to file.
/// Tries to guess the landmark format from the file extension.
void write_landmarks(
        const std::vector<VH>& _landmarks,
        const fs::path& _path,
        const std::vector<LandmarkType>& _landmark_types = std::vector<LandmarkType>(),
        const bool _overwrite = false);

void write_landmarks_pinned(
        const std::vector<VH>& _landmarks,
        const fs::path& _path,
        const std::vector<LandmarkType>& _landmark_types = std::vector<LandmarkType>(),
        const bool _overwrite = false);
void write_landmarks_vts(
        const std::vector<VH>& _landmarks,
        const fs::path& _path,
        const bool _overwrite = false);

/// Read edge indices from file
std::vector<EH> read_paths(
        const std::string& _path,
        const bool _throw_if_not_exist = true);

std::vector<pm::edge_handle> read_paths(
        const std::string& _path,
        const pm::Mesh& _mesh,
        const bool _throw_if_not_exist = true);

/// Write edge indices to file
void write_paths(
        const std::vector<EH>& _edges,
        const std::string& _path,
        const bool _overwrite = false);

void write_paths(
        const std::vector<pm::edge_handle>& _edges,
        const std::string& _path,
        const bool _overwrite = false);

std::shared_ptr<glow::Texture2D> read_texture(
        const fs::path &_file_path);

ExternalProperty<VH, Vec3d> embedding_from_mesh(
        const TriMesh& _mesh);

/// Read vertex positions from mesh with same connectivity
ExternalProperty<VH, Vec3d> read_embedding(
        const TriMesh& _mesh,
        const fs::path& _path);

/// Write vertex positions to mesh with same connectivity
void write_embedding(
        const TriMesh& _mesh,
        const ExternalProperty<VH, Vec3d>& _embedding,
        const fs::path& _path);

/// Convert OpenMesh to libigl matrix representation
void mesh_to_matrix(
        const TriMesh& _mesh,
        MatXd& _V, MatXi& _F);

/// Convert libigl matrix to OpenMesh
void matrix_to_mesh(
        const MatXd& _V,
        const MatXi& _F,
        TriMesh& _mesh,
        const bool _compute_normals = false);

/// Convert parameterization to libigl matrix
MatXd param_to_matrix(
        const Parametrization& _param,
        const TriMesh& _mesh);

/// Convert libigl matrix to parameterization
Parametrization param_from_matrix(
        const MatXd& _P,
        const TriMesh& _mesh);

/// Convert OpenMesh to polymesh
pm::vertex_attribute<Vec3d> to_polymesh(
        const TriMesh& mesh,
        pm::Mesh& m);

pm::vertex_attribute<tg::pos3> to_polymesh_tg(
        const TriMesh& mesh,
        pm::Mesh& m);

MatXd uv_vector_to_matrix(
        const VecXd& _uv_vec);

Parametrization vector_to_param(
        const VecXd& _x,
        const TriMesh _mesh);

VecXd matrix_to_stacked_vector(
        const MatXd& _mat);

/// Read double per mesh element
ExternalProperty<VH, double> read_per_vertex_data(
        const std::string& _path,
        const TriMesh& _mesh,
        const bool _throw_if_not_exist = true);

/// Write double-valued vertex property
void write_vertex_property(
        const OpenMesh::VPropHandleT<double> _ph,
        const std::string& _path,
        const TriMesh& _mesh,
        const bool _overwrite = false);

void write_vertex_to_point_map(
        const VertexToPointMap& _vtpm,
        const std::string& _path,
        const TriMesh& _source,
        const TriMesh& _target,
        const bool _overwrite = false);

VertexToPointMap read_vertex_to_point_map(
        const std::string& _path,
        const TriMesh& _source,
        const TriMesh& _target,
        const bool _throw_if_not_exist = true);

template <typename T>
void write_matrix(
        const MatX<T>& _M,
        const fs::path& _path,
        const bool _overwrite = false,
        const bool _write_header = true);

MatXd read_matrix(
        const fs::path& _path);

template <typename HandleT>
ExternalProperty<HandleT, double> read_property(
        const fs::path& _path,
        const TriMesh& _mesh,
        const bool _throw_if_not_exist);

template <typename HandleT>
void write_property(
        const ExternalProperty<HandleT, double>& _prop,
        const fs::path& _path,
        const bool _overwrite = false);

std::string folder_name_from_time();

std::string pad_integer(
        const int _i,
        const int _pad = 6);

// True iff _flag occours in _argv
bool flag(const char* _flag, int _argc, char** _argv);

}
