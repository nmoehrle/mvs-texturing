#pragma once

#include <string>
#include <vector>
#include <list>

#include "math/vector.h"
#include "math/matrix.h"

#include "mve/view.h"
#include "mve/image_tools.h"
#include "mve/image_io.h"

#include "Tri.h"

#include "Settings.h"

/** Struct containing the quality and mean color of a face within a view. */
struct ProjectedFaceInfo {
    std::uint16_t view_id;
    float quality;
    math::Vec3f mean_color;
};

/** Struct containing the quality of a face within a view. */
struct ReducedProjectedFaceInfo {
    std::uint16_t view_id;
    float quality;
};

inline
ReducedProjectedFaceInfo reduce(ProjectedFaceInfo info) {
    ReducedProjectedFaceInfo reduced = {info.view_id, info.quality};
    return reduced;
}


/**
  * Class representing a view with specialized functions for texturing.
  */
class TextureView {
    private:
        std::size_t id;

        math::Vec3f pos;
        math::Vec3f viewdir;
        math::Matrix3f projection;
        math::Matrix4f world_to_cam;
        int width;
        int height;
        mve::ByteImage::Ptr image;
        mve::ByteImage::Ptr gradient_magnitude;
        std::vector<bool> validity_mask;


    public:
        /** Returns the id of the TexureView which is consistent for every run. */
        std::size_t get_id(void) const;

        /** Returns the 2D pixel coordinates of the given vertex projected into the view. */
        math::Vec2f get_pixel_coords(math::Vec3f const & vertex) const;
        /** Returns the RGB pixel values [0, 1] for the given vertex projected into the view, calculated by linear interpolation. */
        math::Vec3f get_pixel_values(math::Vec3f const & vertex) const;

        /** Returns whether the pixel location is valid in this view.
          * The pixel location is valid if its inside the visible area and,
          * if a validity mask has been generated, all surrounding (integer coordinate) pixels are valid in the validity mask.
          */
        bool valid_pixel(math::Vec2f pixel) const;

        /** TODO */
        bool inside(math::Vec3f const & v1, math::Vec3f const & v2, math::Vec3f const & v3) const;

        /** Returns the RGB pixel values [0, 1] for the give pixel location. */
        math::Vec3f get_pixel_values(math::Vec2f const & pixel) const;

        /** Constructs a TextureView from the give mve::CameraInfo containing the given image. */
        TextureView(std::size_t id, mve::CameraInfo const & camera, mve::ByteImage::Ptr image);

        /** Returns the position. */
        math::Vec3f get_pos(void) const;
        /** Returns the viewing direction. */
        math::Vec3f get_viewing_direction(void) const;
        /** Returns the width of the encapsulated image. */
        int get_width(void) const;
        /** Returns the height of the encapsulated image. */
        int get_height(void) const;
        /** Returns a reference pointer to the encapsulated image. */
        mve::ByteImage::Ptr get_image(void) const;

        /** Exchange encapsulated image. */
        void bind_image(mve::ByteImage::Ptr new_image);

        /** Generates the validity mask. */
        void generate_validity_mask(void);
        /** Generates the gradient magnitude image for the encapsulated image. */
        void generate_gradient_magnitude(void);

        /** Releases the validity mask. */
        void release_validity_mask(void);
        /** Releases the gradient magnitude image. */
        void release_gradient_magnitude(void);
        /** Releases the encapsulated image. */
        void release_image(void);

        /** Erodes the validity mask by one pixel. */
        void erode_validity_mask(void);

        template <DataTerm T> void
        get_face_info(math::Vec3f const & v1, math::Vec3f const & v2, math::Vec3f const & v3, ProjectedFaceInfo * face_info) const;

        void
        export_triangle(math::Vec3f v1, math::Vec3f v2, math::Vec3f v3, std::string const & filename) const;

        void
        export_validity_mask(std::string const & filename) const;
};


inline std::size_t
TextureView::get_id(void) const {
    return id;
}

inline math::Vec3f
TextureView::get_pos(void) const {
    return pos;
}

inline math::Vec3f
TextureView::get_viewing_direction(void) const {
    return viewdir;
}

inline int
TextureView::get_width(void) const {
    return width;
}

inline int
TextureView::get_height(void) const {
    return height;
}

inline mve::ByteImage::Ptr
TextureView::get_image(void) const {
    assert(image != NULL);
    return image;
}

inline bool
TextureView::inside(math::Vec3f const & v1, math::Vec3f const & v2, math::Vec3f const & v3) const {
    math::Vec2f p1 = get_pixel_coords(v1);
    math::Vec2f p2 = get_pixel_coords(v2);
    math::Vec2f p3 = get_pixel_coords(v3);
    return valid_pixel(p1) && valid_pixel(p2) && valid_pixel(p3);
}

inline math::Vec2f
TextureView::get_pixel_coords(math::Vec3f const & vertex) const {
    math::Vec3f pixel = projection * world_to_cam.mult(vertex, 1.0f);
    pixel /= pixel[2];
    return math::Vec2f(pixel[0] - 0.5f, pixel[1] - 0.5f);
}

inline math::Vec3f
TextureView::get_pixel_values(math::Vec3f const & vertex) const {
    math::Vec2f pixel = get_pixel_coords(vertex);
    return get_pixel_values(pixel);
}

inline math::Vec3f
TextureView::get_pixel_values(math::Vec2f const & pixel) const {
    assert(image != NULL);
    math::Vec3uc values;
    image->linear_at(pixel[0], pixel[1], *values);
    return math::Vec3f(values) / 255.0f;
}

inline void
TextureView::bind_image(mve::ByteImage::Ptr new_image) {
    image = new_image;
}

inline void
TextureView::release_validity_mask(void) {
    assert(validity_mask.size() == static_cast<std::size_t>(width * height));
    validity_mask.clear();
    validity_mask.shrink_to_fit();
}

inline void
TextureView::release_gradient_magnitude(void) {
    assert(gradient_magnitude != NULL);
    gradient_magnitude.reset();
}

inline void
TextureView::release_image(void) {
    assert(image != NULL);
    image.reset();
}

template <DataTerm T> void
TextureView::get_face_info(math::Vec3f const & v1, math::Vec3f const & v2,
    math::Vec3f const & v3, ProjectedFaceInfo * face_info) const {

    assert(image != NULL);
    assert(gradient_magnitude != NULL);

    math::Vec2f p1 = get_pixel_coords(v1);
    math::Vec2f p2 = get_pixel_coords(v2);
    math::Vec2f p3 = get_pixel_coords(v3);

    assert(valid_pixel(p1) && valid_pixel(p2) && valid_pixel(p3));

    Tri tri(p1, p2, p3);
    Rect<float> aabb = tri.get_aabb();
    float area = tri.get_area();

    if (area < std::numeric_limits<float>::epsilon()) {
        face_info->quality = 0.0f;
        return;
    }

    math::Vec3d mean_color(0.0);
    math::Accum<math::Vec3f> color_accum(math::Vec3d(0.0));

    /* Only relevant for GMI. */
    double integral = 0.0;

    std::size_t num_samples = 0;
    if (area > 0.5f) {
        /* Sort pixel in ascending order of y */
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

                if (!fast_sampling_possible && !tri.inside(static_cast<float>(x) + 0.5f, static_cast<float>(y) + 0.5f)) continue;

                for (std::size_t i = 0; i < 3; i++){
                     color[i] = static_cast<double>(image->at(x, y, i)) / 255.0;
                }
                color_accum.add(color, 1.0f);

                /* Calculate DataTerms. */
                if (T == GMI) {
                    integral += (static_cast<double>(gradient_magnitude->at(x, y, 0)) / 255.0);
                }
                ++num_samples;
            }
        }
    }

    if (num_samples > 0) {
        mean_color = color_accum.normalized();
        integral = (integral / num_samples) * static_cast<double>(area);
    } else {
        math::Vec3d c1, c2, c3;
        for (std::size_t i = 0; i < 3; ++i) {
             c1[i] = static_cast<double>(image->linear_at(p1[0], p1[1], i)) / 255.0;
             c2[i] = static_cast<double>(image->linear_at(p2[0], p2[1], i)) / 255.0;
             c3[i] = static_cast<double>(image->linear_at(p3[0], p3[1], i)) / 255.0;
        }
        mean_color = ((c1 + c2 + c3) / 3.0);

        if (T == GMI) {
            double gmv1 = static_cast<double>(gradient_magnitude->linear_at(p1[0], p1[1], 0)) / 255.0;
            double gmv2 = static_cast<double>(gradient_magnitude->linear_at(p2[0], p2[1], 0)) / 255.0;
            double gmv3 = static_cast<double>(gradient_magnitude->linear_at(p3[0], p3[1], 0)) / 255.0;
            integral = ((gmv1 + gmv2 + gmv3) / 3.0) * static_cast<double>(area);
        }
    }

    face_info->mean_color = mean_color;
    switch (T) {
        case AREA:
            face_info->quality = area;
        break;
        case GMI:
            face_info->quality = integral;
        break;
    }
}
