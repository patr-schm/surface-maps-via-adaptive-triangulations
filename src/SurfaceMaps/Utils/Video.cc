/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt
 */

#include "Video.hh"
#include <SurfaceMaps/Utils/Out.hh>

namespace SurfaceMaps
{

void render_video(
        const fs::path& _path,
        const int _fps,
        const std::string& _filename)
{
    // ffmpeg -r 30 -f image2 -i '%*.png' -vcodec libx264 -pix_fmt yuv420p -y animation.mp4

    // Convert png files to mp4
    const auto ffmpeg = "ffmpeg -r " + std::to_string(_fps) + " -f image2 -i "
            + (_path / "'%*.png'").string()
            + " -vcodec libx264 -pix_fmt yuv420p -y "
            + (_path / _filename).string();

    ISM_INFO(ffmpeg);
    std::system(ffmpeg.c_str());
}

}
