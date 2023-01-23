/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt
 */

#include "Diffusion.hh"

#include <SurfaceMaps/Utils/IO.hh>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>

namespace SurfaceMaps
{

template <typename T>
void diffuse_pointwise_field(
        const TriMesh& _mesh,
        ExternalProperty<VH, T>& _field,
        const double _t)
{
    MatXd V;
    MatXi F;
    mesh_to_matrix(_mesh, V, F);

    SparseMatrix L;
    igl::cotmatrix(V, F, L);

    SparseMatrix M;
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);

    // Store field as #V-by-k matrix. Rowwise.
    const int k = _field[VH(0)].size();
    MatXd U(_mesh.n_vertices(), k);
    for (auto v : _mesh.vertices())
    {
        for (int row = 0; row < _field[v].rows(); ++row)
        {
            for (int col = 0; col < _field[v].cols(); ++col)
            {
                const int i = row * _field[v].cols() + col;
                U(v.idx(), i) = _field[v](row, col);
            }
        }
    }

    // One step of implicit-Euler integration of heat equation
    // Solve (M - delta * L) U = M * U
    // Based on https://github.com/libigl/libigl/blob/main/tutorial/205_Laplacian/main.cpp
    const auto & S = (M - _t * L);
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double >> solver(S);
    ISM_ASSERT(solver.info() == Eigen::Success);
    U = solver.solve(M * U).eval();

    // Write matrix U back to property
    for (auto v : _mesh.vertices())
    {
        for (int row = 0; row < _field[v].rows(); ++row)
        {
            for (int col = 0; col < _field[v].cols(); ++col)
            {
                const int i = row * _field[v].cols() + col;
                _field[v](row, col) = U(v.idx(), i);
            }
        }
    }
}

// Explicit instantiation
template void diffuse_pointwise_field(const TriMesh&, ExternalProperty<VH, Mat3d>&, const double _t);

void diffuse_pointwise_field(
        const TriMesh& _mesh,
        ExternalProperty<VH, double>& _field,
        const double _t)
{
    ExternalProperty<VH, Mat<1,1,double>> matrix(_mesh);
    for(auto vh : _mesh.vertices())
    {
        matrix[vh](0,0) = _field[vh];
    }
    diffuse_pointwise_field(_mesh, matrix, _t);
    for(auto vh : _mesh.vertices())
    {
         _field[vh] = matrix[vh](0,0);
    }
}

}
