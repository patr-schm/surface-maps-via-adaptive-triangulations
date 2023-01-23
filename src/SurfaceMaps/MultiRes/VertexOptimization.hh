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
#include <SurfaceMaps/MultiRes/MultiResSphereEmbedding.hh>

namespace SurfaceMaps
{

void optimize_vertex(
        MultiResSphereEmbedding& _data,
        const VH _vh,
        const uint _max_iters,
        const MultiResSphereEmbeddingSettings& _settings);

}
