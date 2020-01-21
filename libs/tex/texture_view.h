/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef TEX_TEXTUREVIEW_HEADER
#define TEX_TEXTUREVIEW_HEADER

#include <string>
#include <vector>
#include <list>

#include <math/vector.h>
#include <mve/camera.h>
#include <mve/image.h>
#include <mve/image_io.h>
#include <mve/image_tools.h>

#include "tri.h"
#include "settings.h"

TEX_NAMESPACE_BEGIN

/** Struct containing the quality and mean color of a face within a view. */
struct FaceProjectionInfo {
    std::uint16_t view_id;
    float quality;
    math::Vec3f mean_color;

    bool operator<(FaceProjectionInfo const & other) const {
        return view_id < other.view_id;
    }
};

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
        std::string image_file;
        mve::ImageBase::Ptr image;
        mve::ImageBase::Ptr gradient_magnitude;
        std::vector<bool> validity_mask;
        bool grayscale;

    public:
        /** Returns the id of the TexureView which is consistent for every run. */
        std::size_t get_id(void) const;

        /** Returns the 2D pixel coordinates of the given vertex projected into the view. */
        math::Vec2f get_pixel_coords(math::Vec3f const & vertex) const;
        /** Returns the RGB pixel values [0, 1] for the given vertex projected into the view, calculated by linear interpolation. */
//        math::Vec3f get_pixel_values(math::Vec3f const & vertex) const;

        /** Returns whether the pixel location is valid in this view.
          * The pixel location is valid if its inside the visible area and,
          * if a validity mask has been generated, all surrounding (integer coordinate) pixels are valid in the validity mask.
          */
        bool valid_pixel(math::Vec2f pixel) const;

        /** TODO */
        bool inside(math::Vec3f const & v1, math::Vec3f const & v2, math::Vec3f const & v3) const;

        /** Returns the RGB pixel values [0, 1] for the give pixel location. */
//        math::Vec3f get_pixel_values(math::Vec2f const & pixel) const;

        /** Constructs a TextureView from the give mve::CameraInfo containing the given image. */
        TextureView(std::size_t id, mve::CameraInfo const & camera, std::string const & image_file);

        /** Returns the position. */
        math::Vec3f get_pos(void) const;
        /** Returns the viewing direction. */
        math::Vec3f get_viewing_direction(void) const;
        /** Returns the width of the corresponding image. */
        int get_width(void) const;
        /** Returns the height of the corresponding image. */
        int get_height(void) const;
        /** Returns a reference pointer to the corresponding image. */
        template <typename T>
        typename mve::Image<T>::Ptr
        get_image(void) const;
        mve::ImageBase::Ptr
        get_image(void) const;

        /** Returns a reference pointer to the magnitude gradient. */
        template <typename T>
        typename mve::Image<T>::Ptr
        get_gradient_magnitude(void) const;

        /** Exchange encapsulated image. */
        void bind_image(mve::ImageBase::Ptr new_image);

        /** Loads the corresponding image. */
        void load_image(void);
        /** Generates the validity mask. */
        template <typename T>
        void generate_validity_mask(void);
        /** Generates the gradient magnitude image for the encapsulated image. */
        template <typename T>
        void generate_gradient_magnitude(void);

        /** Releases the validity mask. */
        void release_validity_mask(void);
        /** Releases the gradient magnitude image. */
        void release_gradient_magnitude(void);
        /** Releases the corresponding image. */
        void release_image(void);

        /** Erodes the validity mask by one pixel. */
        void erode_validity_mask(void);

        template <typename T>
        void
        get_face_info(math::Vec3f const & v1, math::Vec3f const & v2, math::Vec3f const & v3,
            FaceProjectionInfo * face_info, Settings const & settings) const;

        void
        export_triangle(math::Vec3f v1, math::Vec3f v2, math::Vec3f v3, std::string const & filename) const;

        void
        export_validity_mask(std::string const & filename) const;

        bool
        is_grayscale();
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

template <typename T>
inline typename mve::Image<T>::Ptr
TextureView::get_image(void) const {
    assert(image != NULL);
    typename mve::Image<T>::Ptr img = std::dynamic_pointer_cast<mve::Image<T>>(image);
    return img;
}

inline mve::ImageBase::Ptr
TextureView::get_image(void) const{
    assert(image != NULL);
    return image;
}

template <typename T>
inline typename mve::Image<T>::Ptr
TextureView::get_gradient_magnitude(void) const {
    if (gradient_magnitude == NULL) return NULL;
    typename mve::Image<T>::Ptr gm = std::dynamic_pointer_cast<mve::Image<T>>(gradient_magnitude);
    return gm;
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

//inline math::Vec3f
//TextureView::get_pixel_values(math::Vec3f const & vertex) const {
//    math::Vec2f pixel = get_pixel_coords(vertex);
//    return get_pixel_values(pixel);
//}

//inline math::Vec3f
//TextureView::get_pixel_values(math::Vec2f const & pixel) const {
//    assert(image != NULL);
//    math::Vec3uc values;
//    image->linear_at(pixel[0], pixel[1], *values);
//    return math::Vec3f(values) / 255.0f;
//}

inline void
TextureView::bind_image(mve::ImageBase::Ptr new_image) {
    image = new_image;
}

inline void
TextureView::release_validity_mask(void) {
    assert(validity_mask.size() == static_cast<std::size_t>(width * height));
    validity_mask = std::vector<bool>();
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

template <typename T>
void
TextureView::generate_validity_mask(void) {
    assert(image != NULL);
    validity_mask.resize(width * height, true);
    mve::ByteImage::Ptr checked = mve::ByteImage::create(width, height, 1);
    typename mve::Image<T>::Ptr image = get_image<T>();

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

template <typename T>
void
TextureView::get_face_info(math::Vec3f const & v1, math::Vec3f const & v2,
    math::Vec3f const & v3, FaceProjectionInfo * face_info, Settings const & settings) const {

    assert(image != NULL);
    assert(settings.data_term != DATA_TERM_GMI || gradient_magnitude != NULL);
    typename mve::Image<T>::Ptr image = get_image<T>();
    typename mve::Image<T>::Ptr gradient_magnitude = get_gradient_magnitude<T>();


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

template <typename T>
void
TextureView::generate_gradient_magnitude(void) {
    assert(image != NULL);
    typename mve::Image<T>::Ptr bw = mve::image::desaturate<T>(get_image<T>(), mve::image::DESATURATE_LUMINANCE);
    gradient_magnitude = mve::image::sobel_edge<T>(bw);
}

inline bool
TextureView::is_grayscale(){
    return grayscale;
}

TEX_NAMESPACE_END

#endif /* TEX_TEXTUREVIEW_HEADER */
