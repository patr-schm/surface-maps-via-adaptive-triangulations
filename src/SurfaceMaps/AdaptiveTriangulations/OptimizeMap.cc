/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Dörte Pieper, Patrick Schmidt
 */

#include <SurfaceMaps/AdaptiveTriangulations/OptimizeMap.hh>

#include <TinyAD/Support/OpenMesh.hh>
#include <TinyAD/ScalarFunction.hh>
#include <TinyAD/Utils/NewtonDirection.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>
#include <TinyAD/Utils/LineSearch.hh>
#include <TinyAD/Utils/Timer.hh>
#include <TinyAD/Utils/Out.hh>

#include <SurfaceMaps/AdaptiveTriangulations/Helpers.hh>
#include <SurfaceMaps/AdaptiveTriangulations/LiftToSurface.hh>
#include <SurfaceMaps/AdaptiveTriangulations/Visualization.hh>
#include <SurfaceMaps/AdaptiveTriangulations/OptimizeMapEnergies.hh>
#include <SurfaceMaps/AdaptiveTriangulations/AssignVerticesToFaces.hh>
#include <SurfaceMaps/AdaptiveTriangulations/AdaptiveTargetEdgeLength.hh>
#include <SurfaceMaps/Utils/IO.hh>

#include <Eigen/IterativeLinearSolvers>

namespace SurfaceMaps
{

OptimizeMapStatus optimize_map(
        MapState& _map_state,
        const int& _iterations,
        std::function<void(const std::string&)> _callback,
        const fs::path& _csv_path,
        const AdaptiveTriangulationsSettings& _settings)
{
    int n_meshes = _map_state.meshes_input.size();
    auto return_value = normal;
    auto num_verts = _map_state.mesh_T.n_vertices();

    // Compute an orthonormal tangent space basis at each vertex.
    std::vector<ExternalProperty<VH, Vec3d>> B1;
    std::vector<ExternalProperty<VH, Vec3d>> B2;
    for (int i = 0; i < n_meshes; i++)
    {
        B1.push_back(ExternalProperty<VH, Vec3d> (_map_state.mesh_T));
        B2.push_back(ExternalProperty<VH, Vec3d> (_map_state.mesh_T));
        compute_local_bases(_map_state.embeddings_T[i], B1[i], B2[i], _map_state.mesh_T);
    }

    // Compute constraint elimination matrix C,
    // which maps from a reduced subspace to the original solution space.
    SparseMatrix C;
    const int n = 2 * n_meshes * _map_state.mesh_T.n_vertices();
    const int m = _settings.hard_landmark_constraints ? n - 2 * n_meshes * _map_state.landmarks_T.size() : n;
    if (_settings.hard_landmark_constraints)
    {
        // Init constrained vertex property
        ExternalProperty<VH, bool> constrained(_map_state.mesh_T, false);
        for (const VH vh_constr : _map_state.landmarks_T)
            constrained[vh_constr] = true;

        // Matrix B maps from m-dim to n-dim space.
        // It is the identity map for all unconstrained vertices.
        std::vector<Triplet> C_triplets;
        C_triplets.reserve(m);
        int C_cols = 0;
        for (const VH vh : _map_state.mesh_T.vertices())
        {
            for (int i = 0; i < n_meshes; i++)
            {
                if (!constrained[vh])
                {
                    C_triplets.emplace_back(2 * vh.idx() + (2 * i * num_verts) + 0, C_cols++, 1.0);
                    C_triplets.emplace_back(2 * vh.idx() + (2 * i * num_verts) + 1, C_cols++, 1.0);
                }
            }
        }
        ISM_ASSERT_EQ(C_cols, m);

        C = SparseMatrix(n, m);
        C.setFromTriplets(C_triplets.cbegin(), C_triplets.cend());
    }
    else
    {
        C = TinyAD::identity<double>(2 * n_meshes * _map_state.mesh_T.n_vertices());
    }
    ISM_ASSERT_EQ(C.rows(), n);
    ISM_ASSERT_EQ(C.cols(), m);

    // Define a retraction operator: A map from a local tangent space to the sphere.
    auto retract = [&] (const auto& v_tang, auto v_idx, int mesh_num)
    {
        // Evaluate target point in 3D ambient space and project to sphere via normalization.
        return (_map_state.embeddings_T[mesh_num][v_idx] + v_tang[0] * B1[mesh_num][v_idx] + v_tang[1] * B2[mesh_num][v_idx]).normalized().eval();
    };

    // Set up function with one 2D tangent vector for each of the n * #v vertices.
    auto func = TinyAD::scalar_function<2>(TinyAD::range(n_meshes  * _map_state.mesh_T.n_vertices()));

    // Add triangle elements for barrier and surface approximation terms.
    func.add_elements<3>(TinyAD::range(n_meshes * _map_state.mesh_T.n_faces()), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element) //for barrier, regularization, surface approx
    {
        // Evaluate element using either double or TinyAD::Double
        using T = TINYAD_SCALAR_TYPE(element);

        // Calculate mesh index of current element
        const int mesh_idx = element.handle / _map_state.mesh_T.n_faces();

        //
        const int triangle_idx = element.handle % _map_state.mesh_T.n_faces();

        const SFH fh = SFH(triangle_idx, &_map_state.mesh_T);

        // Get triangle and vertex handles
        SVH v_a = fh.halfedge().to();
        SVH v_b = fh.halfedge().next().to();
        SVH v_c = fh.halfedge().from();

        // Get 2D vector for mesh_N
        // Get vertex positions, each in their own 2D tangent space.
        // These are (0, 0) when computing derivatives and != (0, 0) in line search.
        Eigen::Vector2<T> a_tang = element.variables(v_a.idx() + (mesh_idx * num_verts));
        Eigen::Vector2<T> b_tang = element.variables(v_b.idx() + (mesh_idx * num_verts));
        Eigen::Vector2<T> c_tang = element.variables(v_c.idx() + (mesh_idx * num_verts));

        // Apply retraction operator: Translate tangent vectors into 3D points on the sphere.
        Eigen::Vector3<T> a_mani_A = retract(a_tang, v_a, mesh_idx);
        Eigen::Vector3<T> b_mani_A = retract(b_tang, v_b, mesh_idx);
        Eigen::Vector3<T> c_mani_A = retract(c_tang, v_c, mesh_idx);

        return eval_singlemesh_energy_triangle_T(
                    a_mani_A, b_mani_A, c_mani_A,
                    mesh_idx,
                    _map_state, _map_state.maps_Tf_inputvs[mesh_idx][fh],
                    _settings);
    });

    // Add triangle pair elements for map distortion.
    if (!_map_state.pairs_map_distortion.empty() && _settings.w_map > 0)
    {
        func.add_elements<6>(TinyAD::range(_map_state.pairs_map_distortion.size()  * _map_state.mesh_T.n_faces()), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element) //for map energy
        {
            // Evaluate element using either double or TinyAD::Double
            using T = TINYAD_SCALAR_TYPE(element);

            // Calculate pair index of current element
            const int pair_idx = element.handle / _map_state.mesh_T.n_faces();

            const int mesh_A_idx = _map_state.pairs_map_distortion[pair_idx].first;
            const int mesh_B_idx = _map_state.pairs_map_distortion[pair_idx].second;
            const int triangle_idx = element.handle % _map_state.mesh_T.n_faces();

            const SFH fh = SFH(triangle_idx, &_map_state.mesh_T);

            // Get triangle and vertex handles
            SVH v_a = fh.halfedge().to();
            SVH v_b = fh.halfedge().next().to();
            SVH v_c = fh.halfedge().from();

            // Get 2D vector for mesh_N
            // Get vertex positions, each in their own 2D tangent space.
            // These are (0, 0) when computing derivatives and != (0, 0) in line search.
            Eigen::Vector2<T> a_tang_A = element.variables(v_a.idx() + (mesh_A_idx * num_verts));
            Eigen::Vector2<T> b_tang_A = element.variables(v_b.idx() + (mesh_A_idx * num_verts));
            Eigen::Vector2<T> c_tang_A = element.variables(v_c.idx() + (mesh_A_idx * num_verts));
            Eigen::Vector2<T> a_tang_B = element.variables(v_a.idx() + (mesh_B_idx * num_verts));
            Eigen::Vector2<T> b_tang_B = element.variables(v_b.idx() + (mesh_B_idx * num_verts));
            Eigen::Vector2<T> c_tang_B = element.variables(v_c.idx() + (mesh_B_idx * num_verts));

            // Apply retraction operator: Translate tangent vectors into 3D points on the sphere.
            Eigen::Vector3<T> a_mani_A = retract(a_tang_A, v_a, mesh_A_idx);
            Eigen::Vector3<T> b_mani_A = retract(b_tang_A, v_b, mesh_A_idx);
            Eigen::Vector3<T> c_mani_A = retract(c_tang_A, v_c, mesh_A_idx);
            Eigen::Vector3<T> a_mani_B = retract(a_tang_B, v_a, mesh_B_idx);
            Eigen::Vector3<T> b_mani_B = retract(b_tang_B, v_b, mesh_B_idx);
            Eigen::Vector3<T> c_mani_B = retract(c_tang_B, v_c, mesh_B_idx);

            return eval_trianglepair_energy_triangle_T(
                        a_mani_A, b_mani_A, c_mani_A,
                        a_mani_B, b_mani_B, c_mani_B,
                        mesh_A_idx, mesh_B_idx,
                        _map_state,
                        _settings);
        });
    }

    // Add triangle pair elements for mesh energy.
    if (_settings.w_mesh > 0.0)
    {
        func.add_elements<6>(TinyAD::range(n_meshes  * _map_state.mesh_T.n_faces()), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
        {
            // Evaluate element using either double or TinyAD::Double
            using T = TINYAD_SCALAR_TYPE(element);

            // Calculate mesh index of current element
            const int mesh_idx = element.handle / _map_state.mesh_T.n_faces();

            //
            const int triangle_idx = element.handle % _map_state.mesh_T.n_faces();

            const SFH fh = SFH(triangle_idx, &_map_state.mesh_T);

            // Get triangle and vertex handles
            SVH v_a = fh.halfedge().to();
            SVH v_b = fh.halfedge().next().to();
            SVH v_c = fh.halfedge().from();

            // Find minimal index for target edge length
            int mesh_tel_idx;
            double min_tel = INFINITY;
            for (int i = 0; i < n_meshes; i++)
            {
                Eigen::Vector2d a_tang_B = element.variables_passive(v_a.idx() + (i * num_verts));
                Eigen::Vector2d b_tang_B = element.variables_passive(v_b.idx() + (i * num_verts));
                Eigen::Vector2d c_tang_B = element.variables_passive(v_c.idx() + (i * num_verts));
                Eigen::Vector3d a_mani_B = retract(a_tang_B, v_a, i);
                Eigen::Vector3d b_mani_B = retract(b_tang_B, v_b, i);
                Eigen::Vector3d c_mani_B = retract(c_tang_B, v_c, i);
                double tel_i = ideal_triangle_edge_length(
                                    a_mani_B, b_mani_B, c_mani_B,
                                    _map_state.meshes_input[i], _map_state.meshes_embeddings_input[i],
                                    _map_state.tels_input[i], _map_state.bsp_embeddings_input[i]);
                if (tel_i < min_tel)
                {
                    mesh_tel_idx = i;
                    min_tel = tel_i;
                }
            }
            ISM_ASSERT_GEQ(mesh_tel_idx, 0);
            ISM_ASSERT_L(mesh_tel_idx, n_meshes);

            // Get 2D vector for mesh_N
            // Get vertex positions, each in their own 2D tangent space.
            // These are (0, 0) when computing derivatives and != (0, 0) in line search.
            Eigen::Vector2<T> a_tang_A = element.variables(v_a.idx() + (mesh_idx * num_verts));
            Eigen::Vector2<T> b_tang_A = element.variables(v_b.idx() + (mesh_idx * num_verts));
            Eigen::Vector2<T> c_tang_A = element.variables(v_c.idx() + (mesh_idx * num_verts));
            Eigen::Vector2<T> a_tang_B = element.variables(v_a.idx() + (mesh_tel_idx * num_verts));
            Eigen::Vector2<T> b_tang_B = element.variables(v_b.idx() + (mesh_tel_idx * num_verts));
            Eigen::Vector2<T> c_tang_B = element.variables(v_c.idx() + (mesh_tel_idx * num_verts));

            // Apply retraction operator: Translate tangent vectors into 3D points on the sphere.
            Eigen::Vector3<T> a_mani_A = retract(a_tang_A, v_a, mesh_idx);
            Eigen::Vector3<T> b_mani_A = retract(b_tang_A, v_b, mesh_idx);
            Eigen::Vector3<T> c_mani_A = retract(c_tang_A, v_c, mesh_idx);
            Eigen::Vector3<T> a_mani_B = retract(a_tang_B, v_a, mesh_tel_idx);
            Eigen::Vector3<T> b_mani_B = retract(b_tang_B, v_b, mesh_tel_idx);
            Eigen::Vector3<T> c_mani_B = retract(c_tang_B, v_c, mesh_tel_idx);

            const T tel = ideal_triangle_edge_length(
                        a_mani_B, b_mani_B, c_mani_B,
                        _map_state.meshes_input[mesh_tel_idx], _map_state.meshes_embeddings_input[mesh_tel_idx],
                        _map_state.tels_input[mesh_tel_idx], _map_state.bsp_embeddings_input[mesh_tel_idx]);

            return eval_mesh_energy_triangle_T(
                        a_mani_A, b_mani_A, c_mani_A,
                        mesh_idx, tel,
                        _map_state,
                        _settings);
        });
    }

    // Add landmark elements
    if (!_settings.hard_landmark_constraints && (_settings.w_landmark_world > 0.0 || _settings.w_landmark_sphere > 0.0))
    {
        func.add_elements<1>(TinyAD::range(n_meshes * _map_state.landmarks_T.size()), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
        {
            // Evaluate element using either double or TinyAD::Double
            using T = TINYAD_SCALAR_TYPE(element);

            // Calculate mesh index of current element
            const int mesh_idx = element.handle / _map_state.landmarks_T.size();
            const int landmark_idx = element.handle % _map_state.landmarks_T.size();

            // Get vertex handles
            const SVH v_T = SVH(_map_state.landmarks_T[landmark_idx].idx(), &_map_state.mesh_T);
            const VH vh_mesh = _map_state.landmarks_input[mesh_idx][landmark_idx];

            // Get 2D vector, first two for A, last two for B
            // Get vertex positions, each in their own 2D tangent space.
            // These are (0, 0) when computing derivatives and != (0, 0) in line search.
            Eigen::Vector2<T> v_tang = element.variables(v_T.idx() + (mesh_idx * num_verts));

            // Apply retraction operator: Translate tangent vectors into 3D points on the sphere.
            Eigen::Vector3<T> v_mani = retract(v_tang, v_T, mesh_idx);

            return eval_landmark_T_energy(v_mani, vh_mesh, mesh_idx, _map_state, _settings);
        });
    }

    // Variable vector: 2D tangent vector per vertex per mesh.
    // Initially, each vertex sits at the tangent-space origin.
    Eigen::VectorXd x = Eigen::VectorXd::Zero(2 * n_meshes * num_verts);

    // Optimize via Projected-Newton
    int iter = 0;
    for (; iter < _iterations; ++iter)
    {
        // Do not store prefactorization because sparsity pattern changes (mesh energy)
        TinyAD::LinearSolver solver;

        TinyAD::Timer timer_iteration("Iteration");

        std::vector<ExternalProperty<VH, Vec3d>> temp_embeddings = _map_state.embeddings_T;

        for (auto v_idx : _map_state.mesh_T.vertices())
        {
            for (int i = 0; i < n_meshes; i++)
                temp_embeddings[i][v_idx] = retract(x.segment(2 * v_idx.idx() + (2 * i * num_verts), 2).eval(), v_idx, i);
        }

        if (_settings.w_approx > 0)
        {
            for (int i = 0; i < n_meshes; i++)
                update_assignment_vertices_to_T_faces(_map_state.maps_Tf_inputvs[i], _map_state.mesh_T, _map_state.embeddings_input[i], temp_embeddings[i]);
        }
        TINYAD_ASSERT(_settings.w_approx <= 0.0 || sanity_check_assignments(_map_state, temp_embeddings));

        TinyAD::Timer timer_derivatives("Derivatives");

        double f; // Objective value
        Eigen::VectorXd g; // Gradient
        Eigen::SparseMatrix<double> H_proj; // Positive-definite Hessian
        func.eval_with_hessian_proj(x, f, g, H_proj, _settings.hessian_projection_eps);
        TINYAD_ASSERT_EQ(g.rows(), n);
        TINYAD_ASSERT_EQ(H_proj.rows(), n);
        TINYAD_ASSERT_EQ(H_proj.cols(), n);

        timer_derivatives.stop();

        TINYAD_DEBUG_OUT("Energy in iteration eval() " << func.eval(x));

        // Compute update direction
        TinyAD::Timer timer_solve("Solve");

        H_proj.makeCompressed(); // Speedup?

        const double nnz_per_row = (double)H_proj.nonZeros() / (double)H_proj.rows();
        const bool use_cg_solver = nnz_per_row > _settings.use_cg_above_nnz_per_row;
        ISM_INFO("Computing Newton step via "
                 << (use_cg_solver ? ("CG tolerance " + std::to_string(_settings.cg_solver_tolerance)) : "direct solve") << ". "
                 << "n = " << x.size() << ", "
                 << "nnz_per_row = " << nnz_per_row);

        Eigen::VectorXd d;
        if (use_cg_solver)
        {
            const Eigen::SparseMatrix<double> H_reduced = C.transpose() * H_proj * C + _settings.w_identity * TinyAD::identity<double>(m);
            Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg_solver;
            cg_solver.setTolerance(_settings.cg_solver_tolerance);
            const Eigen::VectorX<double> d_reduced = cg_solver.compute(H_reduced).solve(-C.transpose() * g);
            TINYAD_ASSERT_FINITE_MAT(d_reduced);
            d = C * d_reduced;
        }
        else
        {
            d = TinyAD::newton_direction_reduced_basis(g, H_proj, C, solver, _settings.w_identity);
        }

        timer_solve.stop();

        // If d is (numerically) not a descent direction, replace it by negative gradient.
        if (d.dot(g) >= 0.0)
        {
            ISM_WARNING("Replacing Newton direction by negative gradient.");
            d = -g;
        }

        // Compute Newton decrement (before limiting step length)
        const double newton_decr = TinyAD::newton_decrement(d, g);
        TINYAD_DEBUG_OUT("newton decr: " << newton_decr);

        // Limit max step size
        double tangent_length_limit = tan(75.0 * M_PI / 180.0); // 75° vertex movement

        if (_map_state.mesh_T.n_vertices() >= 4000)
            tangent_length_limit = tan(10.0 * M_PI / 180.0); // 10° vertex movement

        double max_tangent = 0.0;
        for (int i = 0; i < d.size() / 2; ++i)
            max_tangent = std::max(max_tangent, d.segment(2 * i, 2).norm());
        if (max_tangent > tangent_length_limit)
            d *= tangent_length_limit / max_tangent;

        ISM_DEBUG_VAR(d.norm());

        // Assert descent direction
        ISM_ASSERT_L(d.dot(g), 0.0);

        // Line Search
        TinyAD::Timer timer_line_search("Line search");
        Eigen::VectorXd x_old = x;
        x = TinyAD::line_search(x, d, f, g, [&](const VecXd & x_new)
        {
            for (auto v_idx : _map_state.mesh_T.vertices())
            {
                for (int i = 0; i < n_meshes; i++)
                    temp_embeddings[i][v_idx] = retract(x_new.segment(2 * v_idx.idx() + (2 * i * num_verts), 2).eval(), v_idx, i);
            }
            if (_settings.w_approx > 0)
            {
                for (int i = 0; i < n_meshes; i++)
                    update_assignment_vertices_to_T_faces(_map_state.maps_Tf_inputvs[i], _map_state.mesh_T, _map_state.embeddings_input[i], temp_embeddings[i]);
            }
            TINYAD_ASSERT(_settings.w_approx <= 0.0 || sanity_check_assignments(_map_state, temp_embeddings));
            TINYAD_ASSERT(_settings.w_approx <= 0.0 || sanity_check_all_vertices_assigned(_map_state));

            return func.eval(x_new);
        });

        const double f_new = func.eval(x);

        bool line_search_converged = (x ==  x_old);
        timer_line_search.stop();

        // Re-center local bases
        for (auto v_idx : _map_state.mesh_T.vertices())
        {
            for (int i = 0; i < n_meshes; i++)
                _map_state.embeddings_T[i][v_idx] = retract(x.segment(2 * v_idx.idx() + (2 * i * num_verts), 2).eval(), v_idx, i);
        }

        //TINYAD_ASSERT(_settings.surface_approx_energy_weight <= 0.0 || sanity_check_assignments(_map_state));
        //TINYAD_ASSERT(_settings.surface_approx_energy_weight <= 0.0 || sanity_check_all_vertices_assigned(_map_state));

        for (int i = 0; i < n_meshes; i++)
        {
            compute_local_bases(_map_state.embeddings_T[i], B1[i], B2[i], _map_state.mesh_T);
        }
        x = Eigen::VectorXd::Zero(2 * n_meshes * num_verts);

        timer_iteration.stop();

        // callback
        _callback("after_optimization");

        // Write timings to csv
        if (!_csv_path.empty())
        {
//                "iteration", "iter_type", "derivatives (seconds)", "solve (seconds)", "line search (seconds)", "iteration (seconds)", "objective",
//                "splits (seconds)", "n_splits", "collapses (seconds)", "n_collapses", "flips (seconds)", "n_flips", "n_vertices"
            append_to_csv(_csv_path, iter, "optimize", timer_derivatives.seconds(), timer_solve.seconds(), timer_line_search.seconds(), timer_iteration.seconds(), f, "", "", "", "", "", "", "");
        }

        if (newton_decr < _settings.conv_thresh_newton_decrement)
        {
            ISM_HIGHLIGHT("Newton decrement reached convergence threshold. Stop continuous iterations.");
            return_value = under_threshold;
            iter++;
            break;
        }

//        const double improvement = f - f_new;
//        if (improvement < _settings.conv_thresh_improvement)
//        {
//            ISM_HIGHLIGHT("Improvement reached convergence threshold. Stop continuous iterations.");
//            return_value = under_threshold;
//            iter++;
//            break;
//        }

        if (line_search_converged)
        {
            ISM_HIGHLIGHT("Line search couldn't find improvement. Stop continuous iterations.");
            return_value = line_search_converge;
            iter++;
            if (_settings.w_approx > 0)
            {
                for (int i = 0; i < n_meshes; i++)
                {
                    update_assignment_vertices_to_T_faces(_map_state.maps_Tf_inputvs[i], _map_state.mesh_T, _map_state.embeddings_input[i], _map_state.embeddings_T[i]);
                }
            }
            break;
        }

    }

    TINYAD_ASSERT(_settings.w_approx <= 0.0 || sanity_check_assignments(_map_state));
    TINYAD_ASSERT(_settings.w_approx <= 0.0 || sanity_check_all_vertices_assigned(_map_state));

    TINYAD_DEBUG_OUT("Final energy: " << func.eval(x));

    return return_value;
}

}
