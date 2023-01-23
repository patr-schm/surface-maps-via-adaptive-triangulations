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

namespace SurfaceMaps
{

struct OptimizeSphereEmbeddingSettings
{
    double w_barrier = 1000.0;
    double w_angle = 100.0;
    double w_area = 1.0;
    double max_edge_length_degrees = 140.0;

    double hessian_proj_eps = 1e-9;
    double w_identity = 1e-9;
    double convergence_eps = 1e-3; // Newton decrement threshold
};

bool sphere_embedding_bijective(
        const TriMesh& _mesh,
        const ExternalProperty<VH, Vec3d>& _embedding);

void optimize_sphere_embedding(
        const TriMesh& _mesh,
        ExternalProperty<VH, Vec3d>& _embedding,
        const int _n_iters,
        const OptimizeSphereEmbeddingSettings& _settings,
        std::function<void()> _callback = [] () {});

}
