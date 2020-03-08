/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <math/matrix.h>

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
TextureView::load_image(void) {
    if(image != NULL) return;

    grayscale = false;
    try {
        image = mve::image::load_file(image_file);
    } catch (...) {}

    if (image == NULL){
        try{
            image = mve::image::load_tiff_16_file(image_file);
        } catch (...) {}
    }

    if (image == NULL){
        image = mve::image::load_tiff_float_file(image_file);
    }

    // Assure images have always at least 3 channels
    if (image->channels() == 1){
        if (image->get_type() == mve::IMAGE_TYPE_FLOAT){
            image = mve::image::expand_grayscale<float>(std::dynamic_pointer_cast<mve::FloatImage>(image));
        }else if (image->get_type() == mve::IMAGE_TYPE_UINT16){
            image = mve::image::expand_grayscale<uint16_t>(std::dynamic_pointer_cast<mve::RawImage>(image));
        }else{
            image = mve::image::expand_grayscale<uint8_t>(std::dynamic_pointer_cast<mve::ByteImage>(image));
        }

        grayscale = true;
    }
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

    if (image->get_type() == mve::IMAGE_TYPE_FLOAT){
        mve::image::save_tiff_float_file(mve::image::crop(get_image<float>(), width, height, left, top,
            *math::Vec3us(3.402823466E38, 0, 3.402823466E38)), filename);
    }else if (image->get_type() == mve::IMAGE_TYPE_UINT16){
        mve::image::save_tiff_16_file(mve::image::crop(get_image<uint16_t>(), width, height, left, top,
            *math::Vec3us(65535, 0, 65535)), filename);
    }else{
        mve::image::save_png_file(mve::image::crop(get_image<uint8_t>(), width, height, left, top,
            *math::Vec3uc(255, 0, 255)), filename);
    }
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
