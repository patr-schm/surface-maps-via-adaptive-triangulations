/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt, Janis Born
 */
#pragma once

#include <experimental/filesystem>

namespace SurfaceMaps
{

namespace fs = std::experimental::filesystem;
inline fs::path SOURCE_PATH = fs::path(SOURCE_PATH_STR);
inline fs::path DATA_PATH = fs::path(DATA_PATH_STR);
inline fs::path OUTPUT_PATH = fs::path(OUTPUT_PATH_STR);

}
