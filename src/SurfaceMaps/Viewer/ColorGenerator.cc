/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt, Janis Born
 */

#include "ColorGenerator.hh"

#include <SurfaceMaps/Viewer/Colors.hh>

namespace SurfaceMaps
{

namespace
{

static Color color_cycle[] =
{
    BLUE_100,
    MAGENTA_100,
    GREEN_100,
    RED_100,
    PURPLE_100,
    YELLOW_100,
    PETROL_100,
    ORANGE_100,
    TEAL_100,
    MAY_GREEN_100,
    BORDEAUX_100,
    LILAC_100,

    BLUE_75,
    MAGENTA_75,
    GREEN_75,
    RED_75,
    PURPLE_75,
    YELLOW_75,
    PETROL_75,
    ORANGE_75,
    TEAL_75,
    MAY_GREEN_75,
    BORDEAUX_75,
    LILAC_75,

    BLUE_50,
    MAGENTA_50,
    GREEN_50,
    RED_50,
    PURPLE_50,
    YELLOW_50,
    PETROL_50,
    ORANGE_50,
    TEAL_50,
    MAY_GREEN_50,
    BORDEAUX_50,
    LILAC_50,

    BLUE_25,
    MAGENTA_25,
    GREEN_25,
    RED_25,
    PURPLE_25,
    YELLOW_25,
    PETROL_25,
    ORANGE_25,
    TEAL_25,
    MAY_GREEN_25,
    BORDEAUX_25,
    LILAC_25,
};

template <class T, std::size_t N>
constexpr int array_size(const T (&array)[N]) noexcept
{
    return N;
}

}

Color ColorGenerator::generate_next_color()
{
    Color result = color_cycle[current_index];

    // Advance index
    int n = array_size(color_cycle);
    if (cycle_size > 0)
        n = std::min(n, cycle_size);
    current_index = (current_index + 1) % n;
    return result;
}

std::vector<Color> ColorGenerator::generate_next_colors(int n)
{
    std::vector<Color> colors(n);
    for (int i = 0; i < n; ++i)
        colors[i] = generate_next_color();

    return colors;
}

}
