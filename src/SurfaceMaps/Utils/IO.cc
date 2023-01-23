/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt, Janis Born
 */

#include "IO.hh"

#include <polymesh/Mesh.hh>
#include <glow/objects/Texture2D.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <SurfaceMaps/Utils/Helpers.hh>
#include <SurfaceMaps/Utils/BarycentricPoint.hh>

#include <string>
#include <fstream>

namespace SurfaceMaps
{

/// Create directory (recursive)
void make_directory(
        const std::string& _dir_path)
{
    fs::create_directories(_dir_path);
}

void make_file_directory(
        const std::string& _file_path)
{
    fs::create_directories(fs::path(_file_path).parent_path());
}

std::ifstream open_file_read(
        const std::string& _path,
        const bool _throw_if_not_exist)
{
    std::ifstream file(_path);
    if (!file.good() && _throw_if_not_exist)
        ISM_ERROR_throw("Could not open file " << _path);

    return file;
}

void remove_comment(
        std::string& str,
        std::string comment)
{
    str = str.substr(0, str.find_first_of(comment));
}

void trim(
        std::string& str)
{
    str.erase(0, str.find_first_not_of('\r'));
    str.erase(str.find_last_not_of('\r') + 1 );
    str.erase(0, str.find_first_not_of('\n'));
    str.erase(str.find_last_not_of('\n') + 1 );
    str.erase(0, str.find_first_not_of(' '));
    str.erase(str.find_last_not_of(' ') + 1 );
}

std::vector<std::string> read_line(
        std::ifstream& _file,
        const char _delimiter)
{
    std::vector<std::string> tokens;

    std::string line_string;
    if (!std::getline(_file, line_string))
        return tokens;

    // Remove comments
    remove_comment(line_string, "#");
    trim(line_string);

    if (line_string.empty())
        return tokens;

    // Split string into tokens
    std::istringstream line_stream(line_string);
    std::string token;
    while (std::getline(line_stream, token, _delimiter))
    {
        trim(token);
        if (token.empty())
            continue;

        tokens.push_back(token);
    }

    return tokens;
}

void append_line(
        const fs::path& _file_path,
        const std::string& _line,
        const bool _silent)
{
    std::ofstream file(_file_path, std::ofstream::app);
    ISM_ASSERT(file.good());

    file << _line << std::endl;

    if (!_silent)
        ISM_INFO("Wrote (append) " << _file_path);
}

namespace
{

const auto max_double_precision = std::numeric_limits<long double>::max_digits10 + 2;

bool file_exists(
        const std::string& _file_path)
{
    auto p = fs::path(_file_path);
    return fs::exists(p) && !fs::is_directory(p);
}

bool obj_has_texture_coordinates(
        const fs::path& _path)
{
    ISM_ASSERT_EQ(_path.extension(), ".obj");
    auto file = open_file_read(_path, true);

    std::string line;
    while (std::getline(file, line))
    {
        if (line.find("vt") != std::string::npos)
            return true;
    }

    return false;
}

}

TriMesh read_mesh(
        const fs::path& _path,
        const bool _silent)
{
    if (!file_exists(_path))
        ISM_ERROR_throw("Could not open file " << _path);

    TriMesh mesh;

    // Detect whether the mesh comes with texture coordinates before calling request_...
    // Is this really necessary?
    if (_path.extension() == ".om")
    {
        OpenMesh::IO::Options options;
        options += OpenMesh::IO::Options::FaceTexCoord;
        ISM_ASSERT(OpenMesh::IO::read_mesh(mesh, _path, options));
    }
    else if (_path.extension() == ".obj" && obj_has_texture_coordinates(_path))
    {
        mesh.request_halfedge_texcoords2D();
        OpenMesh::IO::Options options;
        options += OpenMesh::IO::Options::FaceTexCoord;
        ISM_ASSERT(OpenMesh::IO::read_mesh(mesh, _path, options));
    }
    else
    {
        ISM_ASSERT(OpenMesh::IO::read_mesh(mesh, _path));
    }

    if (!_silent)
        ISM_INFO("Read " << _path);

    return mesh;
}

void write_mesh(
        const TriMesh& _mesh,
        const fs::path& _file_path)
{
    make_file_directory(_file_path);

    OpenMesh::IO::Options options;
    if (_mesh.has_vertex_texcoords2D())
        options += OpenMesh::IO::Options::VertexTexCoord;
    if (_mesh.has_halfedge_texcoords2D())
        options += OpenMesh::IO::Options::FaceTexCoord;
    if (_mesh.has_vertex_normals())
        options += OpenMesh::IO::Options::VertexNormal;

    const bool success = OpenMesh::IO::write_mesh(_mesh, _file_path, options, max_double_precision);
    if (success)
        ISM_INFO("Wrote " << _file_path)
    else
        ISM_ERROR_throw("Failed to write " << _file_path);
}

TexCoords tex_coords(
        const TriMesh& _mesh)
{
    ISM_ASSERT(_mesh.has_halfedge_texcoords2D());
    return TexCoords(_mesh, _mesh.halfedges().to_vector([&] (auto h)
    {
        return _mesh.texcoord2D(h);
    }));
}

namespace
{

template <typename T>
std::vector<T> read_element_handles(
        const std::string& _path,
        const bool _throw_if_not_exist)
{
    std::vector<T> handles;
    auto file = open_file_read(_path, _throw_if_not_exist);

    if (!file.good())
        return handles;

    std::string line;
    while (std::getline(file, line))
    {
        // Cut off comments and trim
        remove_comment(line, "#");
        trim(line);

        if (line.empty())
            continue;

        handles.push_back(T(std::stoi(line)));
    }

    ISM_INFO("Read " << _path);

    return handles;
}

template <typename T>
void write_element_handles(
        const std::vector<T>& _handles,
        const std::string& _path,
        const bool _overwrite)
{
    if (file_exists(_path) && !_overwrite)
    {
        ISM_ERROR(_path << " already exists. Not overwriting!");
        return;
    }

    std::ofstream file(_path, std::ofstream::trunc);
    for (const auto& h : _handles)
        file << std::to_string(h.idx()) << std::endl;

    ISM_INFO("Wrote " << _path);
}

template <typename OM_Handle, typename PM_Handle>
std::vector<OM_Handle> openmesh_handles(
        const std::vector<PM_Handle>& _pm)
{
    std::vector<OM_Handle> om;
    for (const auto& h : _pm)
        om.push_back(OM_Handle(h.idx.value));

    return om;
}

template <typename PM_Handle, typename OM_Handle>
std::vector<PM_Handle> polymesh_handles(
        const std::vector<OM_Handle>& _om,
        const pm::Mesh& _mesh)
{
    std::vector<PM_Handle> pm;
    for (const auto& h : _om)
        pm.push_back(_mesh.handle_of(typename PM_Handle::index_t(h.idx())));

    return pm;
}

}

LandmarkFormat guess_landmark_format(
        const fs::path& _path)
{
    std::string ext = _path.extension();
    for (char& c : ext)
        c = std::tolower(c);

    if (ext == ".pinned")
        return LandmarkFormat::Pinned;
    else if (ext == ".vts")
        return LandmarkFormat::VTS;
    else if (ext == ".txt")
        return LandmarkFormat::Pinned;
    else
        return LandmarkFormat::Unknown;
}

std::vector<VH> read_landmarks(
        const fs::path& _path,
        const LandmarkType _type,
        const bool _throw_if_not_exist)
{
    const LandmarkFormat format = guess_landmark_format(_path);
    if (format == LandmarkFormat::Pinned)
        return read_landmarks_pinned(_path, _type, _throw_if_not_exist);
    else if (format == LandmarkFormat::VTS)
        return read_landmarks_vts(_path, _throw_if_not_exist);
    else
        ISM_ERROR_throw("Unrecognized landmark format");
}

std::vector<VH> read_landmarks_pinned(
        const fs::path& _path,
        const LandmarkType _type,
        const bool _throw_if_not_exist)
{
    std::vector<VH> landmarks;
    auto file = open_file_read(_path, _throw_if_not_exist);

    if (!file.good())
        return landmarks;

    while (!file.eof())
    {
        const auto line = read_line(file, ' ');
        if (line.empty())
            continue;

        if (line.size() == 1)
            landmarks.push_back(VH(std::stoi(line[0])));
        else if (line.size() == 2)
        {
            if (line[0] == "li")
            {
                if (_type == LandmarkType::Init)
                    landmarks.push_back(VH(std::stoi(line[1])));
            }
            else if (line[0] == "lp")
            {
                if (_type == LandmarkType::Init ||
                    _type == LandmarkType::PreOpt)
                    landmarks.push_back(VH(std::stoi(line[1])));
            }
            else if (line[0] == "lk")
            {
                landmarks.push_back(VH(std::stoi(line[1])));
            }
            else
                ISM_ERROR_throw("Unsupported landmark format");
        }
        else
            ISM_ERROR_throw("Unsupported landmark format");
    }

    ISM_INFO("Read " << _path);

    return landmarks;
}

std::vector<VH> read_landmarks_vts(
        const fs::path& _path,
        const bool _throw_if_not_exist)
{
    auto f = open_file_read(_path, _throw_if_not_exist);
    std::vector<VH> result;
    while (f.good())
    {
        int id;
        float x, y, z; // Unused
        f >> id >> x >> y >> z;
        if (f.good())
            result.push_back(VH(id));
    }
    return result;
}

void write_landmarks(
        const std::vector<VH>& _landmarks,
        const fs::path& _path,
        const std::vector<LandmarkType>& _landmark_types,
        const bool _overwrite)
{
    const LandmarkFormat format = guess_landmark_format(_path);
    if (format == LandmarkFormat::Pinned)
        write_landmarks_pinned(_landmarks, _path, _landmark_types, _overwrite);
    else if (format == LandmarkFormat::VTS)
        write_landmarks_vts(_landmarks, _path, _overwrite);
    else
        ISM_ERROR_throw("Unrecognized landmark format");
}

void write_landmarks_pinned(
        const std::vector<VH>& _landmarks,
        const fs::path& _path,
        const std::vector<LandmarkType>& _landmark_types,
        const bool _overwrite)
{
    ISM_ASSERT(_landmarks.size() == _landmark_types.size() || _landmark_types.empty());

    if (file_exists(_path) && !_overwrite)
    {
        ISM_ERROR(_path << " already exists. Not overwriting!");
        return;
    }

    std::ofstream file(_path, std::ofstream::trunc);
    for (int i = 0; i < _landmarks.size(); ++i)
    {
        if (_landmark_types.size() == _landmarks.size())
        {
            if (_landmark_types[i] == LandmarkType::Init)
                file << "li ";
            else if (_landmark_types[i] == LandmarkType::PreOpt)
                file << "lp ";
            else if (_landmark_types[i] == LandmarkType::Keep)
                file << "lk ";
            else
                ISM_ERROR_throw("Landmark type not supported.");
        }
        file << std::to_string(_landmarks[i].idx()) << std::endl;
    }

    ISM_INFO("Wrote " << _path);
}

void write_landmarks_vts(
        const std::vector<VH>& _landmarks,
        const fs::path& _path,
        const bool _overwrite)
{
    if (file_exists(_path) && !_overwrite)
    {
        ISM_ERROR(_path << " already exists. Not overwriting!");
        return;
    }
    std::ofstream f(_path);
    for (const auto& vh : _landmarks)
    {
        f << vh.idx(); // ID
        f << " 0.0 0.0 0.0\n"; // Dummy position
    }
}

std::vector<EH> read_paths(
        const std::string& _path,
        const bool _throw_if_not_exist)
{
    return read_element_handles<EH>(_path, _throw_if_not_exist);
}

std::vector<pm::edge_handle> read_paths(
        const std::string& _path,
        const pm::Mesh& _mesh,
        const bool _throw_if_not_exist)
{
    return polymesh_handles<pm::edge_handle>(read_paths(_path, _throw_if_not_exist), _mesh);
}

void write_paths(
        const std::vector<EH>& _edges,
        const std::string& _path,
        const bool _overwrite)
{
    write_element_handles(_edges, _path, _overwrite);
}

void write_paths(
        const std::vector<pm::edge_handle>& _edges,
        const std::string& _path,
        const bool _overwrite)
{
    write_element_handles(openmesh_handles<EH>(_edges), _path, _overwrite);
}

glow::SharedTexture2D read_texture(
        const fs::path &_file_path)
{
    return glow::Texture2D::createFromFile(_file_path, glow::ColorSpace::sRGB);
}

ExternalProperty<VH, Vec3d> embedding_from_mesh(
        const TriMesh& _mesh)
{
    ExternalProperty<VH, Vec3d> embedding(_mesh);
    for (auto v : _mesh.vertices())
        embedding[v] = _mesh.point(v);

    return embedding;
}

ExternalProperty<VH, Vec3d> read_embedding(
        const TriMesh& _mesh,
        const fs::path& _path)
{
    const TriMesh embedded_mesh = read_mesh(_path);
    ISM_ASSERT_EQ(embedded_mesh.n_vertices(), _mesh.n_vertices());

    ExternalProperty<VH, Vec3d> embedding(_mesh);
    for (auto v : _mesh.vertices())
        embedding[v] = embedded_mesh.point(v);

    return embedding;
}

void write_embedding(
        const TriMesh& _mesh,
        const ExternalProperty<VH, Vec3d>& _embedding,
        const fs::path& _path)
{
    ISM_ASSERT(_embedding.size_okay(_mesh));

    TriMesh embedded_mesh = _mesh; // copy
    for (auto v : _mesh.vertices())
        embedded_mesh.point(v) = _embedding[v];

    write_mesh(embedded_mesh, _path);
}

void mesh_to_matrix(
        const TriMesh& _mesh,
        MatXd& _V, MatXi& _F)
{
    _V.resize(_mesh.n_vertices(), 3);
    for (const auto vh : _mesh.vertices())
    {
        const auto& p = _mesh.point(vh);
        for (int col = 0; col < 3; ++col)
            _V(vh.idx(), col) = p[col];
    }

    _F.resize(_mesh.n_faces(), 3);
    for (const auto fh : _mesh.faces())
    {
        const auto heh = _mesh.halfedge_handle(fh);
        _F(fh.idx(), 0) = _mesh.to_vertex_handle(heh).idx();
        _F(fh.idx(), 1) = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(heh)).idx();
        _F(fh.idx(), 2) = _mesh.from_vertex_handle(heh).idx();
    }
}

void matrix_to_mesh(
        const MatXd& _V,
        const MatXi& _F,
        TriMesh& _mesh,
        const bool _compute_normals/* = false*/)
{
    _mesh.clear();
    _mesh.reserve(_V.rows(), _V.rows() + _F.rows(), _F.rows());

    for (int i = 0; i < _V.rows(); ++i)
        _mesh.add_vertex(Vec3d(_V(i, 0), _V(i, 1), _V(i, 2)));
    for (int i = 0; i < _F.rows(); ++i)
        _mesh.add_face(VH(_F(i, 0)), VH(_F(i, 1)), VH(_F(i, 2)));

    if (_compute_normals)
    {
        _mesh.request_face_normals();
        _mesh.update_face_normals();
        _mesh.request_vertex_normals();
        _mesh.update_vertex_normals();
    }
}

MatXd param_to_matrix(
        const Parametrization& _param,
        const TriMesh& _mesh)
{
    MatXd M(_mesh.n_vertices(), 2);
    for (auto vh : _mesh.vertices())
    {
        M(vh.idx(), 0) = _param[vh][0];
        M(vh.idx(), 1) = _param[vh][1];
    }

    return M;
}

/// Convert libigl matrix to parameterization
Parametrization param_from_matrix(
        const MatXd& _P,
        const TriMesh& _mesh)
{
    ISM_ASSERT_EQ(_P.rows(), (int)_mesh.n_vertices());
    ISM_ASSERT_EQ(_P.cols(), 2);

    Parametrization param(_mesh);
    for (auto vh : _mesh.vertices())
        param[vh] = Vec2d(_P(vh.idx(), 0), _P(vh.idx(), 1));

    return param;
}

pm::vertex_attribute<Vec3d> to_polymesh(
        const TriMesh& mesh,
        pm::Mesh& m)
{
    for (auto vh : mesh.vertices())
        m.vertices().add();
    for (auto fh : mesh.faces())
    {
        VH vh_a, vh_b, vh_c;
        handles(mesh, fh, vh_a, vh_b, vh_c);
        m.faces().add(m.handle_of(pm::vertex_index(vh_b.idx())),
                      m.handle_of(pm::vertex_index(vh_c.idx())),
                      m.handle_of(pm::vertex_index(vh_a.idx())));
    }

    auto pos = m.vertices().make_attribute<Vec3d>();
    for (auto vh : mesh.vertices())
        pos[m.handle_of(pm::vertex_index(vh.idx()))] = mesh.point(vh);

    return pos;
}

pm::vertex_attribute<tg::pos3> to_polymesh_tg(
        const TriMesh& mesh,
        pm::Mesh& m)
{
    for (auto vh : mesh.vertices())
        m.vertices().add();
    for (auto fh : mesh.faces())
    {
        VH vh_a, vh_b, vh_c;
        handles(mesh, fh, vh_a, vh_b, vh_c);
        m.faces().add(m.handle_of(pm::vertex_index(vh_b.idx())),
                      m.handle_of(pm::vertex_index(vh_c.idx())),
                      m.handle_of(pm::vertex_index(vh_a.idx())));
    }

    auto pos = m.vertices().make_attribute<tg::pos3>();
    for (auto vh : mesh.vertices())
        pos[m.handle_of(pm::vertex_index(vh.idx()))] = tg::pos3(mesh.point(vh));

    return pos;
}

MatXd uv_vector_to_matrix(
        const VecXd& _uv_vec)
{
    MatXd uv_mat(_uv_vec.rows() / 2, 2);

    for (uint i = 0; i < uv_mat.rows(); ++i)
    {
        uv_mat(i, 0) = _uv_vec[2 * i];
        uv_mat(i, 1) = _uv_vec[2 * i + 1];
    }

    return uv_mat;
}

Parametrization vector_to_param(
        const VecXd& _x,
        const TriMesh _mesh)
{
    Parametrization param(_mesh);

    for (uint i = 0; i < _mesh.n_vertices(); ++i)
    {
        const VH vh(i);
        param[vh] = Vec2d(_x[2 * i], _x[2 * i + 1]);
    }

    return param;
}

VecXd matrix_to_stacked_vector(
        const MatXd &_mat)
{
    VecXd x;
    x.resize(2 * _mat.rows());

    for (uint i = 0; i < _mat.rows(); ++i)
    {
        x[i] = _mat(i, 0);
        x[i + _mat.rows()] = _mat(i, 1);
    }
    return x;
}

/// Read double per mesh element
ExternalProperty<VH, double> read_per_vertex_data(
        const std::string& _path,
        const TriMesh& _mesh,
        const bool _throw_if_not_exist)
{
    return read_property<VH>(_path, _mesh, _throw_if_not_exist);
}

/// Write double-valued vertex property
void write_vertex_property(
        const OpenMesh::VPropHandleT<double> _ph,
        const std::string& _path,
        const TriMesh& _mesh,
        const bool _overwrite)
{
    make_file_directory(_path);

    if (file_exists(_path) && !_overwrite)
    {
        ISM_ERROR(_path << " already exists. Not overwriting!");
        return;
    }

    std::ofstream file(_path);
    ISM_ASSERT(file.is_open());
    for (auto vh : _mesh.vertices())
        file << std::to_string(_mesh.property(_ph, vh)) << std::endl;

    ISM_INFO("Wrote " << _path);
}

void write_vertex_to_point_map(
        const VertexToPointMap& _vtpm,
        const std::string& _path,
        const TriMesh& _source,
        const TriMesh& _target,
        const bool _overwrite)
{
    ISM_ASSERT(_vtpm.size_okay(_source));

    make_file_directory(_path);
    if (file_exists(_path) && !_overwrite)
    {
        ISM_ERROR(_path << " already exists. Not overwriting!");
        return;
    }

    std::ofstream file(_path);
    ISM_ASSERT(file.is_open());

    // Write header: Number of source and target vertices.
    // Just as a very basic way to verify compatibility when loading this map.
    file << _source.n_vertices() << "\n";
    file << _target.n_vertices() << "\n";

    // Write map elements.
    for (const auto& source_vh : _source.vertices())
    {
        const BarycentricPoint& bary = _vtpm[source_vh];
        const HEH target_heh = bary.heh();
        const VH target_vh0 = _target.from_vertex_handle(target_heh);
        const VH target_vh1 = _target.to_vertex_handle(target_heh);

        file << target_vh0.idx() << " ";
        file << target_vh1.idx() << " ";
        file << bary.alpha() << " ";
        file << bary.beta() << "\n";
    }

    ISM_INFO("Wrote " << _path);
}

VertexToPointMap read_vertex_to_point_map(
        const std::string& _path,
        const TriMesh& _source,
        const TriMesh& _target,
        const bool _throw_if_not_exist)
{
    auto file = open_file_read(_path, _throw_if_not_exist);

    int expected_source_vertices = -1;
    int expected_target_vertices = -1;
    file >> expected_source_vertices;
    file >> expected_target_vertices;
    ISM_ASSERT_EQ(expected_source_vertices, _source.n_vertices());
    ISM_ASSERT_EQ(expected_target_vertices, _target.n_vertices());

    VertexToPointMap result(_source);
    for (int source_vi = 0; source_vi < expected_source_vertices; ++source_vi)
    {
        int target_vi0 = -1;
        int target_vi1 = -1;
        double alpha = NAN_DOUBLE;
        double beta = NAN_DOUBLE;

        ISM_ASSERT(file.good());
        file >> target_vi0 >> target_vi1 >> alpha >> beta;

        ISM_ASSERT_GEQ(target_vi0, 0);
        ISM_ASSERT_GEQ(target_vi1, 0);
        ISM_ASSERT_L(target_vi0, expected_target_vertices);
        ISM_ASSERT_L(target_vi1, expected_target_vertices);
        ISM_ASSERT_FINITE(alpha);
        ISM_ASSERT_FINITE(beta);

        const VH target_vh0 = _target.vertex_handle(target_vi0);
        const VH target_vh1 = _target.vertex_handle(target_vi1);
        const HEH target_heh = _target.find_halfedge(target_vh0, target_vh1);
        ISM_ASSERT(target_heh.is_valid());

        const VH source_vh = _source.vertex_handle(source_vi);
        const BarycentricPoint bary(target_heh, alpha, beta, _target);
        result[source_vh] = bary;
    }

    ISM_INFO("Read " << _path);
    return result;
}

MatXd read_matrix(
        const fs::path& _path)
{
    auto file = open_file_read(_path, true);

    int rows, cols;
    file >> rows >> cols;

    MatXd result(rows, cols);
    for (int row = 0; row < rows; ++row)
    {
        for (int col = 0; col < cols; ++col)
        {
            double v;
            file >> v;
            result(row, col) = v;
        }
    }
    return result;
}

template <typename T>
void write_matrix(
        const MatX<T>& _M,
        const fs::path& _path,
        const bool _overwrite,
        const bool _write_header)
{
    if (file_exists(_path) && !_overwrite)
    {
        ISM_ERROR(_path << " already exists. Not overwriting!");
        return;
    }
    make_file_directory(_path);

    std::ofstream file(_path);
    ISM_ASSERT(file.is_open());

    if (_write_header)
        file << _M.rows() << " " << _M.cols() << "\n";
    for (int row = 0; row < _M.rows(); ++row)
    {
        for (int col = 0; col < _M.cols(); ++col)
        {
            file << _M(row, col);
            if (col == _M.cols() - 1)
                file << "\n";
            else
                file << " ";
        }
    }
    ISM_INFO("Wrote " << _path);
}
template void write_matrix(const MatX<int>&, const fs::path&, const bool, const bool);
template void write_matrix(const MatX<double>&, const fs::path&, const bool, const bool);

template <typename HandleT>
ExternalProperty<HandleT, double> read_property(
        const fs::path& _path,
        const TriMesh& _mesh,
        const bool _throw_if_not_exist)
{
    auto file = open_file_read(_path, _throw_if_not_exist);
    if (!file.is_open())
        return ExternalProperty<HandleT, double>();

    ExternalProperty<HandleT, double> data(_mesh);

    std::string line;
    int i = 0;
    while (std::getline(file, line))
    {
        // Cut off comments and trim
        remove_comment(line, "#");
        trim(line);

        if (line.empty())
            continue;

        data[HandleT(i)] = std::stod(line);
        ++i;
    }

    ISM_ASSERT_EQ(i, (int)data.size());
    ISM_INFO("Read " << _path);

    return data;
}
template ExternalProperty<VH, double> read_property(const fs::path&, const TriMesh&, const bool);
template ExternalProperty<EH, double> read_property(const fs::path&, const TriMesh&, const bool);
template ExternalProperty<HEH, double> read_property(const fs::path&, const TriMesh&, const bool);
template ExternalProperty<FH, double> read_property(const fs::path&, const TriMesh&, const bool);

template <typename HandleT>
void write_property(
        const ExternalProperty<HandleT, double>& _prop,
        const fs::path& _path,
        const bool _overwrite)
{
    if (file_exists(_path) && !_overwrite)
    {
        ISM_ERROR(_path << " already exists. Not overwriting!");
        return;
    }
    make_file_directory(_path);

    std::ofstream file(_path);
    file << _prop.as_eigen().format(Eigen::IOFormat(Eigen::FullPrecision));

    ISM_INFO("Wrote " << _path);
}
template void write_property(const ExternalProperty<VH, double>&, const fs::path&, const bool);
template void write_property(const ExternalProperty<EH, double>&, const fs::path&, const bool);
template void write_property(const ExternalProperty<HEH, double>&, const fs::path&, const bool);
template void write_property(const ExternalProperty<FH, double>&, const fs::path&, const bool);

std::string folder_name_from_time()
{
    auto t = std::time(nullptr);

    std::stringstream s;
    s << std::put_time(std::localtime(&t), "%Y_%m_%d_%H_%M_%S");

    return s.str();
}

std::string pad_integer(
        const int _i,
        const int _pad)
{
    std::ostringstream s;
    s << std::setw(_pad) << std::setfill('0') << _i;

    return s.str();
}

bool flag(const char* _flag, int _argc, char** _argv)
{
    for (int i = 1; i < _argc; ++i)
    {
        if (strcmp(_argv[i], _flag) == 0)
            return true;
    }

    return false;
}

}
