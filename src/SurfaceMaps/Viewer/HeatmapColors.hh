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

Color log_color(
        const double _val,
        const double _min,
        const double _max,
        const Color& _min_color,
        const Color& _max_color);

template <typename HandleT>
ExternalProperty<HandleT, Color> log_colors(
        const ExternalProperty<HandleT, double>& _values,
        const double _min,
        const double _max,
        const Color& _min_color,
        const Color& _max_color);

template <typename HandleT>
ExternalProperty<HandleT, Color> linear_colors(
        const ExternalProperty<HandleT, double>& _values,
        const double _min,
        const double _max,
        const Color& _min_color,
        const Color& _max_color);

}
