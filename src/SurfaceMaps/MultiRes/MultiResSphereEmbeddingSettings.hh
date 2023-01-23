/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Joe Jakobi, Patrick Schmidt
 */
#pragma once

#include <SurfaceMaps/Types.hh>
#include <SurfaceMaps/MultiRes/ProgressiveMesh.hh>
#include <SurfaceMaps/AdaptiveTriangulations/OptimizeSphereEmbedding.hh>

namespace SurfaceMaps
{

struct MultiResSphereEmbeddingSettings
{
    MultiResSphereEmbeddingSettings()
    {
        final_opt_settings.w_barrier = 1.0;
        final_opt_settings.w_angle = 100.0;
        final_opt_settings.w_area = 1.0;
        final_opt_settings.convergence_eps = 10.0;

        refinement_opt_settings = final_opt_settings;
    }

    ProgressiveMeshOptions decimation_settings;
    OptimizeSphereEmbeddingSettings refinement_opt_settings;
    OptimizeSphereEmbeddingSettings final_opt_settings;

    bool parallel_optimization = true;
    int phase_one_vertex_thresh = 1000;
    int initialization_max_iters = 10;
    int global_iters_phase_one = 20;
    int global_iters_phase_two = 10;
    int global_iters_final = 50;
    double w_identity = 1e-6;
    double optimization_improvment_thresh = 1e-6;
    double optimization_gradient_thresh = 1e-6;
    bool use_phases = true;
    int iters_exponentional = 200; // Deprecated?
    int iters_minimal = 5; // Deprecated?
    double iters_exponent = 1.5; // Deprecated?

    bool try_hard = true; // Try a bunch of different settings on failure
};

}
