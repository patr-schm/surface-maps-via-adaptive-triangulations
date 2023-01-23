/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt, Janis Born
 */

#include "Geometries.hh"

namespace SurfaceMaps
{

Geometry which_geometry(const TriMesh& _mesh)
{
    if (closed_surface(_mesh))
    {
        if (genus(_mesh) == 0)
            return Spherical;
        else if (genus(_mesh) == 1)
            return Planar;
        else if (genus(_mesh) >= 2)
            return Hyperbolic;
        else
            ISM_ERROR_throw("");
    }
    else
    {
        return Planar;
    }
}

}
