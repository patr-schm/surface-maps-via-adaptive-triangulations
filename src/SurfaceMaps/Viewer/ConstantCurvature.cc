/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt, Janis Born
 */

#include "ConstantCurvature.hh"

#include <SurfaceMaps/Utils/IO.hh>
#include <SurfaceMaps/Utils/Helpers.hh>
#include <SurfaceMaps/Viewer/GlowDraw.hh>
#include <SurfaceMaps/Viewer/MeshView.hh>
#include <SurfaceMaps/Misc/ConstantCurvatureEmbedding.hh>

namespace SurfaceMaps
{

gv::detail::raii_view_closer view_ccm_embedding(
        const TriMesh& _mesh,
        const ExternalProperty<VH, Vec3d>& _embedding,
        const Geometry& _geometry,
        const Color& _color,
        const bool _view_manifold,
        const int _n_samples,
        const DrawStyle& _style)
{
    auto v = gv::view();
//    auto style = default_style();

    ISM_ASSERT(_geometry == Spherical); // TODO: Implement others
    const int n_samples_per_180_deg = _n_samples;

    auto bb = gv::config(tg::aabb3(-1, 1));
    if (_view_manifold)
        view_mesh(read_mesh(DATA_PATH / "meshes/genus0/sphere/sphere.obj", true));

    GlowDraw draw;
    auto draw_arc = [&] (Vec3d from, Vec3d to, Color color, const DrawStyle& style)
    {
        const double dist = acos(dot(from, to));
        const int n_samples = std::max(2, (int)(dist / M_PI * n_samples_per_180_deg));

        for (int i = 0; i < n_samples - 1; ++i)
        {
            const double lambda1 = (double)i / (double)(n_samples - 1);
            const double lambda2 = (double)(i + 1) / (double)(n_samples - 1);
            Vec3d p1 = (1.0 - lambda1) * from + lambda1 * to;
            Vec3d p2 = (1.0 - lambda2) * from + lambda2 * to;
            p1 = p1.normalized();
            p2 = p2.normalized();

            draw.line(p1, p2, color, style);
        }
    };

    for (auto f : _mesh.faces())
    {
        VH vh_a, vh_b, vh_c;
        handles(_mesh, f, vh_a, vh_b, vh_c);

        Color color = _color;
        DrawStyle style = _style;
        if (flipped_or_degenerate(_embedding[vh_a], _embedding[vh_b], _embedding[vh_c], _geometry))
        {
            // Flipped
            color = RED;
            style.width *= 1.5;
        }

        draw_arc(_embedding[vh_a], _embedding[vh_b], color, style);
        draw_arc(_embedding[vh_b], _embedding[vh_c], color, style);
        draw_arc(_embedding[vh_c], _embedding[vh_a], color, style);
    }

    draw.view();

    return v;
}

}
