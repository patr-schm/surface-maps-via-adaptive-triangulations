/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt
 */
#pragma once

#include <SurfaceMaps/Utils/Filesystem.hh>

namespace SurfaceMaps
{

void render_video(
        const fs::path& _path,
        const int _fps,
        const std::string& _filename = "animation.mp4");

}
