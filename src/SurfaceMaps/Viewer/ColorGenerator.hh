/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt, Janis Born
 */

#pragma once

#include <SurfaceMaps/Types.hh>

namespace SurfaceMaps
{

class ColorGenerator
{
public:
    ColorGenerator(int _n_skip = 0, int _cycle_size = -1)
        : current_index(_n_skip),
          cycle_size(_cycle_size) { }

    Color generate_next_color();
    std::vector<Color> generate_next_colors(int n);

private:
    int current_index = 0;
    int cycle_size = -1;
};

}
