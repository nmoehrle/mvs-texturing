#include "TexturePatch.h"

TexturePatch::TexturePatch(int _label, std::vector<std::size_t> const & _faces,
    std::vector<math::Vec2f>  const & _texcoords, mve::ByteImage::Ptr _image)
    : label(_label), faces(_faces.begin(), _faces.end()),
    texcoords(_texcoords.begin(), _texcoords.end()), image(_image) {

    validity_mask = mve::ByteImage::create(get_width(), get_height(), 1);
    validity_mask->fill(255);
    blending_mask = mve::ByteImage::create(get_width(), get_height(), 1);
}

TexturePatch::TexturePatch(TexturePatch const & _texture_patch) {
    label = _texture_patch.get_label();
    faces = std::vector<std::size_t>(_texture_patch.get_faces());
    texcoords = std::vector<math::Vec2f>(_texture_patch.get_texcoords());
    image = mve::ByteImage::create(*(_texture_patch.get_image()));
    validity_mask = mve::ByteImage::create(*(_texture_patch.get_validity_mask()));
    if (_texture_patch.blending_mask != NULL) {
        blending_mask = mve::ByteImage::create(*(_texture_patch.get_blending_mask()));
    }
}

const float sqrt_2 = sqrt(2);

void
TexturePatch::adjust_colors(std::vector<math::Vec3f> const & adjust_values) {
    assert(blending_mask != NULL);

    validity_mask->fill(0);

    mve::FloatImage::Ptr iadjust_values = mve::FloatImage::create(get_width(), get_height(), 3);
    for (std::size_t i = 0; i < texcoords.size(); i += 3) {
        math::Vec2f v1 = texcoords[i];
        math::Vec2f v2 = texcoords[i + 1];
        math::Vec2f v3 = texcoords[i + 2];

        Tri tri(v1, v2, v3);

        float area = tri.get_area();
        if (area < std::numeric_limits<float>::epsilon()) continue;

        Rect<float> aabb = tri.get_aabb();
        int const min_x = static_cast<int>(std::floor(aabb.min_x)) - texture_patch_border;
        int const min_y = static_cast<int>(std::floor(aabb.min_y)) - texture_patch_border;
        int const max_x = static_cast<int>(std::ceil(aabb.max_x)) + texture_patch_border;
        int const max_y = static_cast<int>(std::ceil(aabb.max_y)) + texture_patch_border;
        assert(0 <= min_x && max_x <= get_width());
        assert(0 <= min_y && max_y <= get_height());

        for (int y = min_y; y < max_y; ++y) {
            for (int x = min_x; x < max_x; ++x) {

                math::Vec3f bcoords = tri.get_barycentric_coords(x, y);
                bool inside = bcoords.minimum() >= 0.0f;
                if (inside) {
                    assert(x != 0 && y != 0);
                    for (int c = 0; c < 3; ++c)
                        iadjust_values->at(x, y, c) = math::interpolate(
                            adjust_values[i][c], adjust_values[i + 1][c], adjust_values[i + 2][c],
                            bcoords[0], bcoords[1], bcoords[2]);
                    validity_mask->at(x, y, 0) = 255;
                    blending_mask->at(x, y, 0) = 255;
                } else {

                    if (validity_mask->at(x, y, 0) == 255)
                        continue;

                    /* Check whether the pixels distance from the triangle is more than one pixel. */
                    float ha = 2.0f * -bcoords[0] * area / (v2 - v3).norm();
                    float hb = 2.0f * -bcoords[1] * area / (v1 - v3).norm();
                    float hc = 2.0f * -bcoords[2] * area / (v1 - v2).norm();

                    if (ha > sqrt_2 || hb > sqrt_2 || hc > sqrt_2)
                        continue;

                    for (int c = 0; c < 3; ++c)
                        iadjust_values->at(x, y, c) = math::interpolate(
                            adjust_values[i][c], adjust_values[i + 1][c], adjust_values[i + 2][c],
                            bcoords[0], bcoords[1], bcoords[2]);
                    validity_mask->at(x, y, 0) = 255;
                    blending_mask->at(x, y, 0) = 128;
                }
            }
        }
    }

    mve::FloatImage::Ptr adjusted_image = mve::FloatImage::create(this->get_width(), this->get_height(), 3);
    adjusted_image = mve::image::byte_to_float_image(image);
    for (int i = 0; i < adjusted_image->get_pixel_amount(); ++i)
        if (validity_mask->at(i, 0) != 0){
            for (int c = 0; c < 3; ++c)
                adjusted_image->at(i, c) += iadjust_values->at(i, c);
        } else {
            adjusted_image->at(i, 0) = 1.0f;
            adjusted_image->at(i, 1) = 0.0f;
            adjusted_image->at(i, 2) = 1.0f;
        }

    image = mve::image::float_to_byte_image(adjusted_image);
}

bool TexturePatch::valid_pixel(math::Vec2f pixel) const {
    float x = pixel[0];
    float y = pixel[1];

    float const height = static_cast<float>(get_height());
    float const width = static_cast<float>(get_width());

    bool valid = (0.0f <= x && x < width && 0.0f <= y && y < height);
    if (valid && validity_mask != NULL){
        /* Only pixel which can be correctly interpolated are valid. */
        float cx = std::max(0.0f, std::min(width - 1.0f, x));
        float cy = std::max(0.0f, std::min(height - 1.0f, y));
        int const floor_x = static_cast<int>(cx);
        int const floor_y = static_cast<int>(cy);
        int const floor_xp1 = std::min(floor_x + 1, get_width() - 1);
        int const floor_yp1 = std::min(floor_y + 1, get_height() - 1);

        float const w1 = cx - static_cast<float>(floor_x);
        float const w0 = 1.0f - w1;
        float const w3 = cy - static_cast<float>(floor_y);
        float const w2 = 1.0f - w3;

        valid = (w0 * w2 == 0.0f || validity_mask->at(floor_x, floor_y, 0) == 255) &&
                (w1 * w2 == 0.0f || validity_mask->at(floor_xp1, floor_y, 0) == 255) &&
                (w0 * w3 == 0.0f || validity_mask->at(floor_x, floor_yp1, 0) == 255) &&
                (w1 * w3 == 0.0f || validity_mask->at(floor_xp1, floor_yp1, 0) == 255);
    }

    return valid;
}

bool
TexturePatch::valid_pixel(math::Vec2i pixel) const {
    int const x = pixel[0];
    int const y = pixel[1];

    bool valid = (0 <= x && x < get_width() && 0 <= y && y < get_height());
    if (valid && validity_mask != NULL) {
        valid = validity_mask->at(x, y, 0) == 255;
    }

    return valid;
}

math::Vec3f
TexturePatch::get_pixel_value(math::Vec2f pixel) const {
    assert(valid_pixel(pixel));

    math::Vec3uc color;
    image->linear_at(pixel[0], pixel[1], *color);
    return math::Vec3f(color) / 255.0f;
}

void
TexturePatch::set_pixel_value(math::Vec2i pixel, math::Vec3f color) {
    assert(blending_mask != NULL);
    assert(valid_pixel(pixel));

    for (int c = 0; c < 3; ++c) {
        image->at(pixel[0], pixel[1], c) = color[c] * 255.0f;
    }
    blending_mask->at(pixel[0], pixel[1], 0) = 126;
}

void
TexturePatch::blend(mve::ByteImage::ConstPtr _orig) {
    mve::ByteImage::Ptr orig = mve::ByteImage::create(*_orig);
    poisson_blend(orig, blending_mask, image, 1.0f);
}

void
TexturePatch::prepare_blending_mask(std::size_t strip_width){
    int const width = blending_mask->width();
    int const height = blending_mask->height();

    /* Calculate the set of valid pixels at the border of texture patch. */
    PixelSet valid_border_pixels;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            if (validity_mask->at(x, y, 0) == 255) {
                /* Valid border pixels need no invalid neighbours. */
                if (x == 0 || x == width - 1 || y == 0 || y == height - 1) {
                    valid_border_pixels.insert(std::pair<int, int>(x, y));
                    continue;
                }

                bool at_border = false;
                /* Check the direct neighbourhood of all invalid pixels. */
                for (int j = -1; j <= 1 && !at_border; ++j) {
                    for (int i = -1; i <= 1 && !at_border; ++i) {
                        int nx = x + i;
                        int ny = y + j;
                        /* If the valid pixel has a invalid neighbour: */
                        if (0 <= nx && nx < width &&
                            0 <= ny && ny < height &&
                            validity_mask->at(nx, ny, 0) == 0) {

                            /* Add the pixel to the set of valid border pixels. */
                            valid_border_pixels.insert(std::pair<int, int>(x, y));
                            at_border = true;
                        }
                    }
                }
            }
        }
    }

    mve::ByteImage::Ptr inner_pixel = validity_mask->duplicate();

    /* Iteratively erode all border pixels. */
    for (std::size_t i = 0; i < strip_width; ++i){
        PixelVector new_invalid_pixels(valid_border_pixels.begin(), valid_border_pixels.end());
        PixelVector::iterator it;
        valid_border_pixels.clear();

        /* Mark the new invalid pixels invalid in the validity mask. */
        for (it = new_invalid_pixels.begin(); it != new_invalid_pixels.end(); ++it) {
             int x = it->first;
             int y = it->second;

             inner_pixel->at(x, y, 0) = 0;
        }

        /* Calculate the set of valid pixels at the border of the valid area. */
        for (it = new_invalid_pixels.begin(); it != new_invalid_pixels.end(); ++it) {
             int x = it->first;
             int y = it->second;

             for (int j = -1; j <= 1; j++){
                 for (int i = -1; i <= 1; i++){
                     int nx = x + i;
                     int ny = y + j;
                     if (0 <= nx && nx < width &&
                         0 <= ny && ny < height &&
                         inner_pixel->at(nx, ny, 0) == 255){

                         valid_border_pixels.insert(std::pair<int, int>(nx, ny));
                     }
                 }
             }
        }
    }

    /* Mark all remaining pixels invalid in the blending_mask. */
    for (int i = 0; i < inner_pixel->get_pixel_amount(); ++i) {
        if (inner_pixel->at(i) == 255) blending_mask->at(i) = 0;
    }

    /* Mark all border pixels. */
    PixelSet::iterator it;
    for (it = valid_border_pixels.begin(); it != valid_border_pixels.end(); ++it) {
         int x = it->first;
         int y = it->second;

         blending_mask->at(x, y, 0) = 128;
    }
}

void
TexturePatch::erode_validity_mask(void) {

    int const width = validity_mask->width();
    int const height = validity_mask->height();

    mve::ByteImage::Ptr eroded_validity_mask = validity_mask->duplicate();

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            if (x == 0 || x == width - 1 || y == 0 || y == height - 1)
                eroded_validity_mask->at(x, y, 0) = 0;

            bool invalid = validity_mask->at(x, y, 0) == 0;
            for (int j = -1; j <= 1 && !invalid; ++j) {
                for (int i = -1; i <= 1 && !invalid; ++i) {
                    int nx = x + i;
                    int ny = y + j;
                    if (0 <= nx && nx < width &&
                        0 <= ny && ny < height &&
                        validity_mask->at(nx, ny, 0) == 0) {
                        eroded_validity_mask->at(x, y, 0) = 0;
                        invalid = true;
                    }
                }
            }
        }
    }

    validity_mask.swap(eroded_validity_mask);
}
