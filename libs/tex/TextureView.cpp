#include "TextureView.h"

TextureView::TextureView(std::size_t id, mve::CameraInfo const & camera,
    mve::ByteImage::Ptr image)
    : id(id), image(image) {

    width = image->width();
    height = image->height();

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

            bool invalid = !validity_mask[x + y * width];
            for (int j = -1; j <= 1 && !invalid; ++j) {
                for (int i = -1; i <= 1 && !invalid; ++i) {
                    int nx = x + i;
                    int ny = y + j;
                    if (0 <= nx && nx < width &&
                        0 <= ny && ny < height &&
                        !validity_mask[nx + ny * width]) {
                        eroded_validity_mask[x + y * width] = false;
                        invalid = true;
                    }
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
