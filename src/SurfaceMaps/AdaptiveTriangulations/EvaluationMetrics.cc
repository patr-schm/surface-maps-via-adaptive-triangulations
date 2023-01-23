/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Dörte Pieper, Patrick Schmidt
 */

#include <SurfaceMaps/AdaptiveTriangulations/EvaluationMetrics.hh>

#include <SurfaceMaps/AdaptiveTriangulations/Helpers.hh>
#include <SurfaceMaps/AdaptiveTriangulations/LiftToSurface.hh>
#include <SurfaceMaps/Utils/Helpers.hh>
#include <SurfaceMaps/Utils/IO.hh>
#include <SurfaceMaps/Utils/BSPTree.hh>
#include <SurfaceMaps/Utils/MeshNormalization.hh>

#include <Eigen/SVD>

#include <igl/hausdorff.h>
#include <igl/point_mesh_squared_distance.h>

namespace SurfaceMaps
{

TriMesh map_vertices_to_target(
        const MapState& _map_state,
        const int _source_mesh_idx,
        const int _target_mesh_idx)
{
    TriMesh mesh_on_target = _map_state.meshes_input[_source_mesh_idx];
    TriMesh mesh_embedding_T_source = embedding_to_mesh(_map_state.mesh_T, _map_state.embeddings_T[_source_mesh_idx]);
    BSPTree bsp_T_source(mesh_embedding_T_source);

    for (auto vh : _map_state.meshes_input[_source_mesh_idx].vertices())
    {
        // Compute barycentric coordinates of points of source mesh in T
        SFH fh_T;
        double alpha_T, beta_T, gamma_T;
        bsp_tree_barys_face(_map_state.embeddings_input[_source_mesh_idx][vh], mesh_embedding_T_source, bsp_T_source, alpha_T, beta_T, gamma_T, fh_T);
        BarycentricPoint bary_T(fh_T, alpha_T, beta_T, _map_state.mesh_T);

        // Map points of source mesh onto target mesh sphere
        const Vec3d p_emb_on_ref = bary_T.interpolate(_map_state.embeddings_T[_target_mesh_idx], _map_state.mesh_T).normalized();

        mesh_on_target.point(vh) = lift_vertex_to_surface(p_emb_on_ref, _map_state.meshes_input[_target_mesh_idx], _map_state.meshes_embeddings_input[_target_mesh_idx], _map_state.bsp_embeddings_input[_target_mesh_idx]);
    }

    return mesh_on_target;
}

double avg_min_angle_for_mesh(const TriMesh& _mesh)
{
    double sum_A = 0.0;
    for (auto fh : _mesh.faces())
    {
        auto heh0 = fh.halfedge();
        auto heh1 = heh0.next();
        auto heh2 = heh1.next();
        auto alpha = _mesh.calc_sector_angle(heh0);
        auto beta = _mesh.calc_sector_angle(heh1);
        auto gamma = _mesh.calc_sector_angle(heh2);
        sum_A += std::min({alpha, beta, gamma});
    }
    //average
    sum_A = sum_A / _mesh.n_faces();
    //to degree
    return sum_A * 180.0 / EIGEN_PI;
}

//hausdorff distance points mesh_1 -> surface mesh_2
double compute_hausdorff_distance(
        const TriMesh& _mesh_1,
        const TriMesh& _mesh_2,
        const int& _samples_per_face)
{
    double distance = 0.0;
    BSPTree bsp_2(_mesh_2);

    MatXd V2;
    MatXi F2;
    mesh_to_matrix(_mesh_2, V2, F2);
    //ISM_DEBUG_OUT(V2.rows() << " " << V2.cols() << " " << _mesh_2.n_vertices());
    //ISM_DEBUG_OUT(_mesh_1.n_vertices() << " " << _mesh_1.n_edges() << " " << _mesh_1.n_faces());

    VecXd sqrD;
    VecXi I;
    MatXd C;

    const int row_num = _mesh_1.n_vertices() + _mesh_1.n_edges() +_mesh_1.n_faces() * _samples_per_face;
    Eigen::MatrixXd P(row_num, 3);
    //ISM_DEBUG_OUT(P.rows() << " " << P.cols());
    int i = 0;

    //vertices
    for (auto vh1 : _mesh_1.vertices())
    {
        P.row(i) = _mesh_1.point(vh1).transpose();
        i++;
    }
    //edges (middle point of edges)
    for (auto eh1 : _mesh_1.edges())
    {
        auto heh = eh1.h0();
        auto vec = _mesh_1.calc_edge_vector(heh);
        P.row(i) = (_mesh_1.point(heh.from()) + 0.5 * vec).transpose();
        i++;
    }
    //faces
    for (auto fh1 : _mesh_1.faces())
    {
        Vec3d p0, p1, p2;
        points(_mesh_1, fh1, p0, p1, p2);
        auto vec1 = p1-p0;
        auto vec2 = p2-p0;
        for (int j =0; j<_samples_per_face; j++){
            auto rand1 = (double)rand() / (double)RAND_MAX;
            auto rand2 = (double)rand() / (double)RAND_MAX;
            if (rand1 + rand2 > 1.0)
            {
                rand1 = 1.0 - rand1;
                rand2 = 1.0 - rand1;
            }
            ISM_ASSERT(rand1 <= 1.0 && rand2 <= 1.0);
            P.row(i) = (p0 + (vec1 * rand1 + vec2 * rand2)).transpose();
            i++;
        }


    }

    igl::point_mesh_squared_distance(P, V2, F2, sqrD, I, C);
    distance = sqrD.maxCoeff();

    return sqrt(distance);
}

double hausdorff_distance_for_trimeshes(const TriMesh& _mesh_orig,
                                        const TriMesh& _mesh_approx)
{
    MatXd V_X, V_T_X;
    MatXi F_X, F_T_X;
    mesh_to_matrix(_mesh_orig, V_X, F_X);
    mesh_to_matrix(_mesh_approx, V_T_X, F_T_X);

    //compute two sided Hausdorff distance between mesh_T_A and _map_state.mesh_A
    double distance;
    igl::hausdorff(V_X, F_X, V_T_X, F_T_X, distance);

    return distance;
}

void mapping_distortion_normalized(
        const TriMesh& _mesh_TA,
        const TriMesh& _mesh_TB,
        ExternalProperty<FH, double>& _areas_TA,
        ExternalProperty<FH, double>& _areas_TB,
        ExternalProperty<FH, double>& _sing_vals_min,
        ExternalProperty<FH, double>& _sing_vals_max)
{
    ISM_ASSERT_EQ(_mesh_TA.n_faces(), _mesh_TB.n_faces());

    _areas_TA.init(_mesh_TA);
    _areas_TB.init(_mesh_TB);
    _sing_vals_min.init(_mesh_TA);
    _sing_vals_max.init(_mesh_TA);

    TriMesh mesh_A_normalized = _mesh_TA;
    TriMesh mesh_B_normalized = _mesh_TB;

    normalize_mesh(mesh_A_normalized);
    normalize_mesh(mesh_B_normalized);

    for (auto fh : _mesh_TA.faces())
    {
        SVH vh_a, vh_b, vh_c;
        handles(_mesh_TA, fh, vh_a, vh_b, vh_c);
        //local coordinate system per face for lifted vertices
        Eigen::Vector2<double> a_lifted_local_A, b_lifted_local_A, c_lifted_local_A, a_lifted_local_B, b_lifted_local_B, c_lifted_local_B;
        to_local_coordinates(mesh_A_normalized.point(vh_a), mesh_A_normalized.point(vh_b), mesh_A_normalized.point(vh_c), a_lifted_local_A, b_lifted_local_A, c_lifted_local_A);
        to_local_coordinates(mesh_B_normalized.point(vh_a), mesh_B_normalized.point(vh_b), mesh_B_normalized.point(vh_c), a_lifted_local_B, b_lifted_local_B, c_lifted_local_B);

        Eigen::Matrix2<double> M_A;
        M_A << b_lifted_local_A - a_lifted_local_A, c_lifted_local_A - a_lifted_local_A;
        Eigen::Matrix2<double> M_B;
        M_B << b_lifted_local_B - a_lifted_local_B, c_lifted_local_B - a_lifted_local_B;

        //Jacobian
        Eigen::Matrix2<double> J = M_B * M_A.inverse();
        Eigen::JacobiSVD<Eigen::Matrix2<double>> svd(J);
        //ISM_DEBUG_OUT("singular values " << svd.singularValues());

        _areas_TA[fh] = mesh_A_normalized.calc_face_area(fh);
        _areas_TB[fh] = mesh_B_normalized.calc_face_area(fh);
        _sing_vals_min[fh] = svd.singularValues()[1];
        _sing_vals_max[fh] = svd.singularValues()[0];
        ISM_ASSERT_LEQ(_sing_vals_min[fh], _sing_vals_max[fh]);
    }
}

void write_mapping_distortion_normalized(
        const TriMesh& _mesh_TA,
        const TriMesh& _mesh_TB,
        const fs::path& _csv_path)
{
    ISM_ASSERT_EQ(_mesh_TA.n_faces(), _mesh_TB.n_faces());

    ExternalProperty<FH, double> areas_TA;
    ExternalProperty<FH, double> areas_TB;
    ExternalProperty<FH, double> sing_vals_min;
    ExternalProperty<FH, double> sing_vals_max;
    mapping_distortion_normalized(_mesh_TA, _mesh_TB, areas_TA, areas_TB, sing_vals_min, sing_vals_max);

    fs::remove_all(_csv_path);
    append_to_csv(_csv_path, "area_TA", "area_TB", "sing_val_min", "sing_val_max");
    for (auto f : _mesh_TA.faces())
        append_to_csv_silent(_csv_path, areas_TA[f], areas_TB[f], sing_vals_min[f], sing_vals_max[f]);
}

void metric_mapping_distortion(const TriMesh& _mesh_1,
                               const TriMesh& _mesh_2,
                               double & _distortion)
{
    double sum_distortion = 0.0;
    for (auto fh : _mesh_1.faces())
    {
        SVH vh_a, vh_b, vh_c;
        handles(_mesh_1, fh, vh_a, vh_b, vh_c);
        //local coordinate system per face for lifted vertices
        Eigen::Vector2<double> a_lifted_local_A, b_lifted_local_A, c_lifted_local_A, a_lifted_local_B, b_lifted_local_B, c_lifted_local_B;
        to_local_coordinates(_mesh_1.point(vh_a), _mesh_1.point(vh_b), _mesh_1.point(vh_c), a_lifted_local_A, b_lifted_local_A, c_lifted_local_A);
        to_local_coordinates(_mesh_2.point(vh_a), _mesh_2.point(vh_b), _mesh_2.point(vh_c), a_lifted_local_B, b_lifted_local_B, c_lifted_local_B);

        Eigen::Matrix2<double> M_A;
        M_A << b_lifted_local_A - a_lifted_local_A, c_lifted_local_A - a_lifted_local_A;
        Eigen::Matrix2<double> M_B;
        M_B << b_lifted_local_B - a_lifted_local_B, c_lifted_local_B - a_lifted_local_B;

        //Jacobian
        Eigen::Matrix2<double> J = M_B * M_A.inverse();

        Eigen::JacobiSVD<Eigen::Matrix2<double>> svd(J);
        //ISM_DEBUG_OUT("singular values " << svd.singularValues());
        double max_distortion = svd.singularValues()[0];
        double min_distortion = svd.singularValues()[1];
        //ISM_DEBUG_OUT("max " << max_distortion << " 1/min " << 1.0/min_distortion);

        //WARNING sometimes result is nan or inf, especially for overlay meshes
        // Prevent nan values in result
        if (isnan(min_distortion) || isnan(max_distortion))
        {
            ISM_WARNING("A distortion value is nan!");
            min_distortion = 1.0;
            max_distortion = 1.0;
        }
        if (isinf(min_distortion) || isinf(max_distortion))
        {
            ISM_WARNING("A distortion value is inf!");
            min_distortion = 1.0;
            max_distortion = 1.0;
        }
        if (iszero(min_distortion)){
            ISM_WARNING("min distortion zero!");
            min_distortion = 1.0;
            max_distortion = 1.0;
        }
        if (min_distortion < 0.0001){
            ISM_WARNING("min distortion small! " << min_distortion);
            min_distortion = 1.0;
        }


        //sum_distortion += (max_distortion + 1.0/min_distortion)/2;
        sum_distortion += fmax((max_distortion), 1.0/(min_distortion)); //TODO seems to be same numbers for table 1 in yang2020
    }
    _distortion = sum_distortion / _mesh_1.n_faces();

}

void evaluate_meshes(const TriMesh& _mesh_A,
                     const TriMesh& _mesh_A_approx,
                     const TriMesh& _mesh_B,
                     const TriMesh& _mesh_B_approx,
                     const double _runtime_seconds,
                     const fs::path& _out_path,
                     const bool& _silent)
{
    double angle_A, angle_B, distortion, distance_A, distance_B;
    angle_A = avg_min_angle_for_mesh(_mesh_A_approx);
    angle_B = avg_min_angle_for_mesh(_mesh_B_approx);
    distance_A = hausdorff_distance_for_trimeshes(_mesh_A, _mesh_A_approx);
    distance_B = hausdorff_distance_for_trimeshes(_mesh_B, _mesh_B_approx);

    //metric_average_min_angles(_map_state, angle_A, angle_B);
    metric_mapping_distortion(_mesh_A_approx, _mesh_B_approx, distortion);
    //metric_hausdorff_distances(_map_state, distance_A, distance_B);

    double dist_test_A = compute_hausdorff_distance(_mesh_A, _mesh_A_approx, 10);
    double dist_test_B = compute_hausdorff_distance(_mesh_B, _mesh_B_approx, 10);
    double dist_test_back_A = compute_hausdorff_distance(_mesh_A_approx, _mesh_A, 10);
    double dist_test_back_B = compute_hausdorff_distance(_mesh_B_approx, _mesh_B, 10);

    double percent_A, percent_B;
    //diagonal length
    double diag_A = length_diagonal_of_bounding_box(_mesh_A);
    double diag_B = length_diagonal_of_bounding_box(_mesh_B);
    ISM_DEBUG_OUT("diagonal A: " << diag_A << " B: " << diag_B);
    percent_A = (distance_A / diag_A) * 100.0;
    percent_B = (distance_B / diag_B) * 100.0;

    if (!_silent)
    {
        ISM_INFO("Number of Vertices: " << _mesh_A_approx.n_vertices() << " B: " << _mesh_B_approx.n_vertices());
        ISM_INFO("Number of Faces: " << _mesh_A_approx.n_faces() << " B: " << _mesh_B_approx.n_faces());
        ISM_INFO("Algorithm took: " <<  _runtime_seconds << " seconds");
        ISM_INFO("Distortion: " << distortion);
        ISM_INFO("Avg smallest angle for A: " << angle_A << " for B: " << angle_B);
        ISM_INFO("Hausdorff distance for A: " << distance_A << " for B: " << distance_B);
        ISM_DEBUG_OUT("other Hausdorff distance for A: " << dist_test_A << " for B: " << dist_test_B);
        ISM_DEBUG_OUT("other Hausdorff distance for A: " << dist_test_back_A << " for B: " << dist_test_back_B);
        ISM_INFO("Hausdorff percentage for A: " << percent_A << "%, for B: " << percent_B <<"%");
    }

    if (_out_path != "")
    {
        fs::remove(_out_path);
        make_file_directory(_out_path);

        //std::string test = _mesh_A_approx.n_vertices() + ", ";
        append_line(_out_path, "Number of Vertices: " + std::to_string(_mesh_A_approx.n_vertices()));
        append_line(_out_path, "Number of Faces: " + std::to_string(_mesh_A_approx.n_faces()));
        append_line(_out_path, "Algorithm took: " + std::to_string(_runtime_seconds) + " seconds");
        append_line(_out_path, "Distortion: " + std::to_string(distortion));
        append_line(_out_path, "Avg smallest angle for A: " + std::to_string(angle_A));
        append_line(_out_path, "Avg smallest angle for B: " + std::to_string(angle_B));
        append_line(_out_path, "Hausdorff percentage for A " + std::to_string(percent_A) + "%");
        append_line(_out_path, "Hausdorff percentage for B " + std::to_string(percent_B) + "%");
    }
}


}
