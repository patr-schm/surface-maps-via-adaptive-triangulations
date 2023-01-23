/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt, Janis Born
 */

#include "Init.hh"

#include <ExactPredicates.h>

namespace SurfaceMaps
{

void init_lib_surface_maps()
{
    exactinit();
}

}
