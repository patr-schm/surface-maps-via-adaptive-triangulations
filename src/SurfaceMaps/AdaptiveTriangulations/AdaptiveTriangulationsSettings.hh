/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Dörte Pieper, Patrick Schmidt
 */
#pragma once

#include <SurfaceMaps/Types.hh>
#include <SurfaceMaps/AdaptiveTriangulations/Helpers.hh>

namespace SurfaceMaps
{

struct AdaptiveTriangulationsSettings
{
    // Schedule
    int max_iterations = 200;
    int optimize_per_remeshing_iters = 1;

    // Weights
    double w_barrier = 1e-4;
    double w_reg = 0.0;
    double w_mesh = 1.0;
    double w_map = 1.0;
    double w_landmark_sphere = 0.0;
    double w_landmark_world = 0.0;
    double w_approx = 1.0;
    double w_bound = 0.0; // Surface approximation bound weight
    double surface_approx_bound = INFINITY;
    bool hard_landmark_constraints = false;
    double w_angle = 0.0;

    // Target edge length
    double approx_error = 0.0025;
    double min_edge_length = 0.0001;
    double max_edge_length = 0.1;

    // Solver settings
    double hessian_projection_eps = 1e-9;
    double w_identity = 1e-9;
    int use_cg_above_nnz_per_row = 50;
    double cg_solver_tolerance = 1e-9;

    // Convergence criteria
    double conv_thresh_newton_decrement = 0.01;
    double landmark_threshold = -INF_DOUBLE; // Stop if remeshing converged and landmark distance in world space is under thresh

    // Remeshing
    bool allow_splits = true;
    bool allow_collapses = true;
    bool allow_flips = true;

    void set_adaptive_tel_params(double approx_e, double min_el, double max_el)
    {
        approx_error = approx_e;
        min_edge_length = min_el;
        max_edge_length = max_el;
    }
};

}
