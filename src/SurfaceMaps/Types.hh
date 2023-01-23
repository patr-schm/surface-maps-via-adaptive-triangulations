/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt, Janis Born
 */
#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Geometry/EigenVectorT.hh>

#include <SurfaceMaps/Utils/ExternalProperty.hh>
#include <SurfaceMaps/Utils/Out.hh>

namespace SurfaceMaps
{

// Scalar constants
constexpr double NAN_DOUBLE = std::numeric_limits<double>::quiet_NaN();
constexpr double INF_DOUBLE = std::numeric_limits<double>::infinity();
constexpr int DEFAULT_MAX_ITERS = -1;

// Complex numbers
using Complex = std::complex<double>;

// Fixed-size vectors
template <int d, typename T> using Vec = Eigen::Matrix<T, d, 1>;
template <typename T> using Vec2 = Vec<2, T>;
template <typename T> using Vec3 = Vec<3, T>;
template <typename T> using Vec4 = Vec<4, T>;
template <typename T> using Vec6 = Vec<6, T>;
template <typename T> using Vec9 = Vec<9, T>;
using Vec2d = Vec2<double>;
using Vec2i = Vec2<int>;
using Vec3d = Vec3<double>;
using Vec3i = Vec3<int>;
using Vec4f = Vec4<float>;
using Vec4d = Vec4<double>;
using Vec6d = Vec6<double>;
using Vec9d = Vec9<double>;

// Fixed-size matrices
template <int rows, int cols, typename T> using Mat = Eigen::Matrix<T, rows, cols>;
template <typename T> using Mat2 = Mat<2, 2, T>;
template <typename T> using Mat3 = Mat<3, 3, T>;
template <typename T> using Mat6 = Mat<6, 6, T>;
template <typename T> using Mat9 = Mat<9, 9, T>;
using Mat2d = Mat2<double>;
using Mat3d = Mat3<double>;
using Mat6d = Mat6<double>;
using Mat9d = Mat9<double>;

// Dynamic-size vectors
template <typename T> using VecX = Eigen::Matrix<T, Eigen::Dynamic, 1>;
using VecXd = Eigen::VectorXd;
using VecXi = Eigen::VectorXi;

// Dynamic-size matrices
template <typename T> using MatX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
using MatXd = MatX<double>;
using MatXi = MatX<int>;

// Sparse matrices
using Triplet = Eigen::Triplet<double>;
using SparseMatrix = Eigen::SparseMatrix<double>;

// Diagonal matrices
using DiagonalMatrix = Eigen::DiagonalMatrix<double, Eigen::Dynamic>;

// Custom vectors
using Color = Vec4f;

// OpenMesh traits
struct ISMOpenMeshTraits : public OpenMesh::DefaultTraits
{
  typedef Vec3d Point;
  typedef Vec3d Normal;
  typedef double TexCoord1D;
  typedef Vec2d TexCoord2D;
  typedef Vec3d TexCoord3D;
  typedef Color Color;
};

// Mesh types
using TriMesh = OpenMesh::TriMesh_ArrayKernelT<ISMOpenMeshTraits>;
using PolyMesh = OpenMesh::PolyMesh_ArrayKernelT<ISMOpenMeshTraits>;

// Handle types
using VH = OpenMesh::VertexHandle;
using EH = OpenMesh::EdgeHandle;
using FH = OpenMesh::FaceHandle;
using HEH = OpenMesh::HalfedgeHandle;
using SVH = OpenMesh::SmartVertexHandle;
using SHEH = OpenMesh::SmartHalfedgeHandle;
using SEH = OpenMesh::SmartEdgeHandle;
using SFH = OpenMesh::SmartFaceHandle;

// Custom properties
using Parametrization = ExternalProperty<VH, Vec2d>;
using Metric = ExternalProperty<EH, double>;
using TexCoords = ExternalProperty<HEH, Vec2d>;

struct BarycentricPoint;
using VertexToPointMap = ExternalProperty<VH, BarycentricPoint>;

}
