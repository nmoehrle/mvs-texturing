/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <list>

#include <math/matrix.h>
#include <mve/image_io.h>
#include <mve/image_tools.h>

#include "texture_view.h"

TEX_NAMESPACE_BEGIN

TextureView::TextureView(std::size_t id, mve::CameraInfo const & camera,
    std::string const & image_file)
    : id(id), image_file(image_file) {

    mve::image::ImageHeaders header;
    try {
         header = mve::image::load_file_headers(image_file);
    } catch (util::Exception e) {
        std::cerr << "Could not load image header of " << image_file << std::endl;
        std::cerr << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    width = header.width;
    height = header.height;

    camera.fill_calibration(*projection, width, height);
    camera.fill_camera_pos(*pos);
    camera.fill_viewing_direction(*viewdir);
    camera.fill_world_to_cam(*world_to_cam);
}

void
TextureView::generate_validity_mask(void) {
    assert(image != NULL);
    validity_mask.resize(width * height, true);
    mve::ByteImage::Ptr checked = mve::ByteImage::create(width, height, 1);

    std::list<math::Vec2i> queue;

    /* Start from the corners. */
    queue.push_back(math::Vec2i(0,0));
    checked->at(0, 0, 0) = 255;
    queue.push_back(math::Vec2i(0, height - 1));
    checked->at(0, height - 1, 0) = 255;
    queue.push_back(math::Vec2i(width - 1, 0));
    checked->at(width - 1, 0, 0) = 255;
    queue.push_back(math::Vec2i(width - 1, height - 1));
    checked->at(width - 1, height - 1, 0) = 255;

    while (!queue.empty()) {
        math::Vec2i pixel = queue.front();
        queue.pop_front();

        int const x = pixel[0];
        int const y = pixel[1];

        int sum = 0;
        for (int c = 0; c < image->channels(); ++c) {
            sum += image->at(x, y, c);
        }

        if (sum == 0) {
            validity_mask[x + y * width] = false;

            std::vector<math::Vec2i> neighbours;
            neighbours.push_back(math::Vec2i(x + 1, y));
            neighbours.push_back(math::Vec2i(x, y + 1));
            neighbours.push_back(math::Vec2i(x - 1, y));
            neighbours.push_back(math::Vec2i(x, y - 1));

            for (std::size_t i = 0; i < neighbours.size(); ++i) {
                math::Vec2i npixel = neighbours[i];
                int const nx = npixel[0];
                int const ny = npixel[1];
                if (0 <= nx && nx < width && 0 <= ny && ny < height) {
                    if (checked->at(nx, ny, 0) == 0) {
                        queue.push_front(npixel);
                        checked->at(nx, ny, 0) = 255;
                    }
                }
            }
        }
    }
}

void
TextureView::load_image(void) {
    if(image != NULL) return;
    image = mve::image::load_file(image_file);
}

void
TextureView::generate_gradient_magnitude(void) {
    assert(image != NULL);
    mve::ByteImage::Ptr bw = mve::image::desaturate<std::uint8_t>(image, mve::image::DESATURATE_LUMINANCE);
    gradient_magnitude = mve::image::sobel_edge<std::uint8_t>(bw);
}

void
TextureView::erode_validity_mask(void) {
    std::vector<bool> eroded_validity_mask(validity_mask);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            if (x == 0 || x == width - 1 || y == 0 || y == height - 1) {
                validity_mask[x + y * width] = false;
                continue;
            }

            if (validity_mask[x + y * width]) continue;
            for (int j = -1; j <= 1; ++j) {
                for (int i = -1; i <= 1; ++i) {
                    int const nx = x + i;
                    int const ny = y + j;
                    eroded_validity_mask[nx + ny * width] = false;
                }
            }
        }
    }

    validity_mask.swap(eroded_validity_mask);
}

void
TextureView::get_face_info(math::Vec3f const & v1, math::Vec3f const & v2,
    math::Vec3f const & v3, FaceProjectionInfo * face_info, Settings const & settings) const {

    assert(image != NULL);
    assert(settings.data_term != DATA_TERM_GMI || gradient_magnitude != NULL);

    math::Vec2f p1 = get_pixel_coords(v1);
    math::Vec2f p2 = get_pixel_coords(v2);
    math::Vec2f p3 = get_pixel_coords(v3);

    assert(valid_pixel(p1) && valid_pixel(p2) && valid_pixel(p3));

    Tri tri(p1, p2, p3);
    float area = tri.get_area();

    if (area < std::numeric_limits<float>::epsilon()) {
        face_info->quality = 0.0f;
        return;
    }

    std::size_t num_samples = 0;
    math::Vec3d colors(0.0);
    double gmi = 0.0;

    bool sampling_necessary = settings.data_term != DATA_TERM_AREA || settings.outlier_removal != OUTLIER_REMOVAL_NONE;

    if (sampling_necessary && area > 0.5f) {
        /* Sort pixels in ascending order of y */
        while (true)
            if(p1[1] <= p2[1])
                if(p2[1] <= p3[1]) break;
                else std::swap(p2, p3);
            else std::swap(p1, p2);

        /* Calculate line equations. */
        float const m1 = (p1[1] - p3[1]) / (p1[0] - p3[0]);
        float const b1 = p1[1] - m1 * p1[0];

        /* area != 0.0f => m1 != 0.0f. */
        float const m2 = (p1[1] - p2[1]) / (p1[0] - p2[0]);
        float const b2 = p1[1] - m2 * p1[0];

        float const m3 = (p2[1] - p3[1]) / (p2[0] - p3[0]);
        float const b3 = p2[1] - m3 * p2[0];

        bool fast_sampling_possible = std::isfinite(m1) && m2 != 0.0f && std::isfinite(m2) && m3 != 0.0f && std::isfinite(m3);

        Rect<float> aabb = tri.get_aabb();
        for (int y = std::floor(aabb.min_y); y < std::ceil(aabb.max_y); ++y) {
            float min_x = aabb.min_x - 0.5f;
            float max_x = aabb.max_x + 0.5f;

            if (fast_sampling_possible) {
                float const cy = static_cast<float>(y) + 0.5f;

                min_x = (cy - b1) / m1;
                if (cy <= p2[1]) max_x = (cy - b2) / m2;
                else max_x = (cy - b3) / m3;

                if (min_x >= max_x) std::swap(min_x, max_x);

                if (min_x < aabb.min_x || min_x > aabb.max_x) continue;
                if (max_x < aabb.min_x || max_x > aabb.max_x) continue;
            }

            for (int x = std::floor(min_x + 0.5f); x < std::ceil(max_x - 0.5f); ++x) {
                math::Vec3d color;

                const float cx = static_cast<float>(x) + 0.5f;
                const float cy = static_cast<float>(y) + 0.5f;
                if (!fast_sampling_possible && !tri.inside(cx, cy)) continue;

                if (settings.outlier_removal != OUTLIER_REMOVAL_NONE) {
                    for (std::size_t i = 0; i < 3; i++){
                         color[i] = static_cast<double>(image->at(x, y, i)) / 255.0;
                    }
                    colors += color;
                }

                if (settings.data_term == DATA_TERM_GMI) {
                    gmi += static_cast<double>(gradient_magnitude->at(x, y, 0)) / 255.0;
                }
                ++num_samples;
            }
        }
    }

    if (settings.data_term == DATA_TERM_GMI) {
        if (num_samples > 0) {
            gmi = (gmi / num_samples) * area;
        } else {
            double gmv1 = static_cast<double>(gradient_magnitude->linear_at(p1[0], p1[1], 0)) / 255.0;
            double gmv2 = static_cast<double>(gradient_magnitude->linear_at(p2[0], p2[1], 0)) / 255.0;
            double gmv3 = static_cast<double>(gradient_magnitude->linear_at(p3[0], p3[1], 0)) / 255.0;
            gmi = ((gmv1 + gmv2 + gmv3) / 3.0) * area;
        }
    }

    if (settings.outlier_removal != OUTLIER_REMOVAL_NONE) {
        if (num_samples > 0) {
            face_info->mean_color = colors / num_samples;
        } else {
            math::Vec3d c1, c2, c3;
            for (std::size_t i = 0; i < 3; ++i) {
                 c1[i] = static_cast<double>(image->linear_at(p1[0], p1[1], i)) / 255.0;
                 c2[i] = static_cast<double>(image->linear_at(p2[0], p2[1], i)) / 255.0;
                 c3[i] = static_cast<double>(image->linear_at(p3[0], p3[1], i)) / 255.0;
            }
            face_info->mean_color = ((c1 + c2 + c3) / 3.0);
        }
    }

    switch (settings.data_term) {
        case DATA_TERM_AREA: face_info->quality = area; break;
        case DATA_TERM_GMI:  face_info->quality = gmi; break;
    }
}

bool
TextureView::valid_pixel(math::Vec2f pixel) const {
    float const x = pixel[0];
    float const y = pixel[1];

    /* The center of a pixel is in the middle. */
    bool valid = (x >= 0.0f && x < static_cast<float>(width - 1)
        && y >= 0.0f && y < static_cast<float>(height - 1));

    if (valid && validity_mask.size() == static_cast<std::size_t>(width * height)) {
        /* Only pixel which can be correctly interpolated are valid. */
        float cx = std::max(0.0f, std::min(static_cast<float>(width - 1), x));
        float cy = std::max(0.0f, std::min(static_cast<float>(height - 1), y));
        int const floor_x = static_cast<int>(cx);
        int const floor_y = static_cast<int>(cy);
        int const floor_xp1 = std::min(floor_x + 1, width - 1);
        int const floor_yp1 = std::min(floor_y + 1, height - 1);

        /* We screw up if weights would be zero
         * e.g. we lose valid pixel in the border of images... */

        valid = validity_mask[floor_x + floor_y * width] &&
                validity_mask[floor_x + floor_yp1 * width] &&
                validity_mask[floor_xp1 + floor_y * width] &&
                validity_mask[floor_xp1 + floor_yp1 * width];
    }

    return valid;
}

void
TextureView::export_triangle(math::Vec3f v1, math::Vec3f v2, math::Vec3f v3,
    std::string const & filename) const {
    assert(image != NULL);
    math::Vec2f p1 = get_pixel_coords(v1);
    math::Vec2f p2 = get_pixel_coords(v2);
    math::Vec2f p3 = get_pixel_coords(v3);

    assert(valid_pixel(p1) && valid_pixel(p2) && valid_pixel(p3));

    Tri tri(p1, p2, p3);

    Rect<float> aabb = tri.get_aabb();
    const int width = ceil(aabb.width());
    const int height = ceil(aabb.height());
    const int left = floor(aabb.min_x);
    const int top = floor(aabb.max_y);

    assert(width > 0 && height > 0);
    mve::image::save_png_file(mve::image::crop(image, width, height, left, top,
        *math::Vec3uc(255, 0, 255)), filename);
}

void
TextureView::export_validity_mask(std::string const & filename) const {
    assert(validity_mask.size() == static_cast<std::size_t>(width * height));
    mve::ByteImage::Ptr img = mve::ByteImage::create(width, height, 1);
    for (std::size_t i = 0; i < validity_mask.size(); ++i) {
        img->at(static_cast<int>(i), 0) = validity_mask[i] ? 255 : 0;
    }
    mve::image::save_png_file(img, filename);
}

TEX_NAMESPACE_END
