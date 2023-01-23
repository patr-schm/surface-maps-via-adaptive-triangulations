/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt, Janis Born
 */

#include "HeatmapColors.hh"

namespace SurfaceMaps
{

Color log_color(
        const double _val,
        const double _min,
        const double _max,
        const Color& _min_color,
        const Color& _max_color)
{
    ISM_ASSERT_G(_val, 0.0);

    double lambda = (std::log(_val) - std::log(_min)) / (std::log(_max) - std::log(_min));
    lambda = std::min(std::max(lambda, 0.0), 1.0);

    return (1.0 - lambda) * _min_color + lambda * _max_color;
}

template <typename HandleT>
ExternalProperty<HandleT, Color> log_colors(
        const ExternalProperty<HandleT, double>& _values,
        const double _min,
        const double _max,
        const Color& _min_color,
        const Color& _max_color)
{
    ExternalProperty<HandleT, Color> colors(_values.size());
    for (int i = 0; i < (int)_values.container_.size(); ++i)
        colors.container_[i] = log_color(_values.container_[i], _min, _max, _min_color, _max_color);

    return colors;
}

// Explicit instantiation
template ExternalProperty<VH, Color> log_colors(const ExternalProperty<VH, double>&, const double, const double, const Color&, const Color&);
template ExternalProperty<FH, Color> log_colors(const ExternalProperty<FH, double>&, const double, const double, const Color&, const Color&);

template <typename HandleT>
ExternalProperty<HandleT, Color> linear_colors(
        const ExternalProperty<HandleT, double>& _values,
        const double _min,
        const double _max,
        const Color& _min_color,
        const Color& _max_color)
{
    ExternalProperty<HandleT, Color> colors(_values.size());
    for (int i = 0; i < (int)_values.container_.size(); ++i)
    {
        double lambda = (_values.container_[i] - _min) / (_max - _min);
        lambda = std::min(std::max(0.0, lambda), 1.0); // clamp to [0..1]
        colors.container_[i] = (1.0 - lambda) * _min_color + lambda * _max_color;
    }

    return colors;
}

// Explicit instantiation
template ExternalProperty<VH, Color> linear_colors(const ExternalProperty<VH, double>&, const double, const double, const Color&, const Color&);
template ExternalProperty<FH, Color> linear_colors(const ExternalProperty<FH, double>&, const double, const double, const Color&, const Color&);
}
