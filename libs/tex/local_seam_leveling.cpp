/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <math/accum.h>

#include "progress_counter.h"
#include "texturing.h"
#include "seam_leveling.h"

TEX_NAMESPACE_BEGIN

#define STRIP_SIZE 20

math::Vec3f
mean_color_of_edge_point(std::vector<EdgeProjectionInfo> const & edge_projection_infos,
    std::vector<TexturePatch::Ptr> const & texture_patches, float t) {

    assert(0.0f <= t && t <= 1.0f);
    math::Accum<math::Vec3f> color_accum(math::Vec3f(0.0f));

    for (EdgeProjectionInfo const & edge_projection_info : edge_projection_infos) {
        TexturePatch::Ptr texture_patch = texture_patches[edge_projection_info.texture_patch_id];
        if (texture_patch->get_label() == 0) continue;
        math::Vec2f pixel = edge_projection_info.p1 * t + (1.0f - t) * edge_projection_info.p2;
        math::Vec3f color = texture_patch->get_pixel_value(pixel);
        if (std::isnan(color[0])) continue;
        color_accum.add(color, 1.0f);
    }

    if (color_accum.w > 0){
        math::Vec3f mean_color = color_accum.normalized();
        return mean_color;
    }else{
        return math::Vec3f(1.0f);
    }
}

void
draw_line(math::Vec2f p1, math::Vec2f p2,
    std::vector<math::Vec3f> const & edge_color, TexturePatch::Ptr texture_patch) {
    /* http://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm */

    int x0 = std::floor(p1[0] + 0.5f);
    int y0 = std::floor(p1[1] + 0.5f);
    int const x1 = std::floor(p2[0] + 0.5f);
    int const y1 = std::floor(p2[1] + 0.5f);

    float tdx = static_cast<float>(x1 - x0);
    float tdy = static_cast<float>(y1 - y0);
    float length = std::sqrt(tdx * tdx + tdy * tdy);

    int const dx = std::abs(x1 - x0);
    int const dy = std::abs(y1 - y0);
    int const sx = x0 < x1 ? 1 : -1;
    int const sy = y0 < y1 ? 1 : -1;
    int err = dx - dy;

    int x = x0;
    int y = y0;
    while (true) {
        math::Vec2i pixel(x, y);

        tdx = static_cast<float>(x1 - x);
        tdy = static_cast<float>(y1 - y);

        /* If the length is zero we sample the midpoint of the projected edge. */
        float t = (length != 0.0f) ? std::sqrt(tdx * tdx + tdy * tdy) / length : 0.5f;

        math::Vec3f color;
        if (t < 1.0f && edge_color.size() > 1) {
            std::size_t idx = std::floor(t * (edge_color.size() - 1));
            color = (1.0f - t) * edge_color[idx] + t * edge_color[idx + 1];
        } else {
            color = edge_color.back();
        }

        texture_patch->set_pixel_value(pixel, color);
        if (x == x1 && y == y1)
            break;

        int const e2 = 2 * err;
        if (e2 > -dy) {
            err -= dy;
            x += sx;
        }
        if (e2 < dx) {
            err += dx;
            y += sy;
        }
    }
}

struct Pixel {
    math::Vec2i pos;
    math::Vec3f const * color;
};

struct Line {
    math::Vec2i from;
    math::Vec2i to;
    std::vector<math::Vec3f> const * color;
};

void
local_seam_leveling(UniGraph const & graph, mve::TriangleMesh::ConstPtr mesh,
    VertexProjectionInfos const & vertex_projection_infos,
    std::vector<TexturePatch::Ptr> * texture_patches) {

    std::size_t const num_vertices = vertex_projection_infos.size();
    std::vector<math::Vec3f> vertex_colors(num_vertices);
    std::vector<std::vector<math::Vec3f> > edge_colors;

    std::vector<std::vector<EdgeProjectionInfo> > edge_projection_infos;
    {
        std::vector<MeshEdge> seam_edges;
        find_seam_edges(graph, mesh, &seam_edges);
        edge_colors.resize(seam_edges.size());
        edge_projection_infos.resize(seam_edges.size());
        for (std::size_t i = 0; i < seam_edges.size(); ++i) {
            MeshEdge const & seam_edge = seam_edges[i];
            find_mesh_edge_projections(vertex_projection_infos, seam_edge,
                &edge_projection_infos[i]);
        }
    }

    std::vector<std::vector<Line> > lines(texture_patches->size());
    std::vector<std::vector<Pixel> > pixels(texture_patches->size());
    /* Sample edge colors. */
    for (std::size_t i = 0; i < edge_projection_infos.size(); ++i) {
        /* Determine sampling (ensure at least two samples per edge). */
        float max_length = 1;
        for (EdgeProjectionInfo const & edge_projection_info : edge_projection_infos[i]) {
            float length = (edge_projection_info.p1 - edge_projection_info.p2).norm();
            max_length = std::max(max_length, length);
        }

        std::vector<math::Vec3f> & edge_color = edge_colors[i];
        edge_color.resize(std::ceil(max_length * 2.0f));
        for (std::size_t j = 0; j < edge_color.size(); ++j) {
            float t = static_cast<float>(j) / (edge_color.size() - 1);
            edge_color[j] = mean_color_of_edge_point(edge_projection_infos[i], *texture_patches, t);
        }

        for (EdgeProjectionInfo const & edge_projection_info : edge_projection_infos[i]) {
            Line line;
            line.from = edge_projection_info.p1 + math::Vec2f(0.5f, 0.5f);
            line.to = edge_projection_info.p2 + math::Vec2f(0.5f, 0.5f);
            line.color = &edge_colors[i];
            lines[edge_projection_info.texture_patch_id].push_back(line);
        }
    }

    /* Sample vertex colors. */
    for (std::size_t i = 0; i < vertex_colors.size(); ++i) {
        std::vector<VertexProjectionInfo> const & projection_infos = vertex_projection_infos[i];
        if (projection_infos.size() <= 1) continue;

        math::Accum<math::Vec3f> color_accum(math::Vec3f(0.0f));
        for (VertexProjectionInfo const & projection_info : projection_infos) {
            TexturePatch::Ptr texture_patch = texture_patches->at(projection_info.texture_patch_id);
            if (texture_patch->get_label() == 0) continue;
            math::Vec3f color = texture_patch->get_pixel_value(projection_info.projection);
            color_accum.add(color, 1.0f);
        }
	if (color_accum.w == 0.0f) continue;

        vertex_colors[i] = color_accum.normalized();

        for (VertexProjectionInfo const & projection_info : projection_infos) {
            Pixel pixel;
            pixel.pos = math::Vec2i(projection_info.projection + math::Vec2f(0.5f, 0.5f));
            pixel.color = &vertex_colors[i];
            pixels[projection_info.texture_patch_id].push_back(pixel);
        }
    }

    ProgressCounter texture_patch_counter("\tBlending texture patches", texture_patches->size());
    #pragma omp parallel for schedule(dynamic)
    for (std::size_t i = 0; i < texture_patches->size(); ++i) {
        TexturePatch::Ptr texture_patch = texture_patches->at(i);
        mve::FloatImage::Ptr image = texture_patch->get_image()->duplicate();

        /* Apply colors. */
        for (Pixel const & pixel : pixels[i]) {
            texture_patch->set_pixel_value(pixel.pos, *pixel.color);
        }

        for (Line const & line : lines[i]) {
            draw_line(line.from, line.to, *line.color, texture_patch);
        }

        texture_patch_counter.progress<SIMPLE>();

        /* Only alter a small strip of texture patches originating from input images. */
        if (texture_patch->get_label() != 0) {
            texture_patch->prepare_blending_mask(STRIP_SIZE);
        }

        texture_patch->blend(image);
        texture_patch->release_blending_mask();
        texture_patch_counter.inc();
    }
}

TEX_NAMESPACE_END
