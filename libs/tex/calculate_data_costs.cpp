/*
 * Copyright (C) 2015, Nils Moehrle, Michael Waechter
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <numeric>

#include <mve/image_color.h>
#include <coldet.h>
#include <Eigen/Core>
#include <Eigen/LU>

#include "util.h"
#include "histogram.h"
#include "texturing.h"
#include "sparse_table.h"
#include "progress_counter.h"

TEX_NAMESPACE_BEGIN

/**
 * Dampens the quality of all views in which the face's projection
 * has a much different color than in the majority of views.
 * Returns whether the outlier removal was successfull.
 *
 * @param infos contains information about one face seen from several views
 * @param settings runtime configuration.
 */
bool
photometric_outlier_detection(std::vector<ProjectedFaceInfo> * infos, Settings const & settings) {
    if (infos->size() == 0) return true;

    /* Configuration variables. */

    double const gauss_rejection_threshold = 6e-3;

    /* If all covariances drop below this we stop outlier detection. */
    double const minimal_covariance = 5e-4;

    int const outlier_detection_iterations = 10;
    int const minimal_num_inliers = 4;

    float outlier_removal_factor = std::numeric_limits<float>::signaling_NaN();
    switch (settings.outlier_removal) {
        case NONE: return true;
        case GAUSS_CLAMPING:
            outlier_removal_factor = 1.0f;
        break;
        case GAUSS_DAMPING:
            outlier_removal_factor = 0.2f;
        break;
    }

    Eigen::MatrixX3d inliers(infos->size(), 3);
    std::vector<std::uint32_t> is_inlier(infos->size(), 1);
    for (std::size_t row = 0; row < infos->size(); ++row) {
        inliers.row(row) = mve_to_eigen(infos->at(row).mean_color).cast<double>();
    }

    Eigen::RowVector3d var_mean;
    Eigen::Matrix3d covariance;
    Eigen::Matrix3d covariance_inv;

    for (int i = 0; i < outlier_detection_iterations; ++i) {

        if (inliers.rows() < minimal_num_inliers) {
            return false;
        }

        /* Calculate the inliers' mean color and color covariance. */
        var_mean = inliers.colwise().mean();
        Eigen::MatrixX3d centered = inliers.rowwise() - var_mean;
        covariance = (centered.adjoint() * centered) / double(inliers.rows() - 1);

        /* If all covariances are very small we stop outlier detection
         * and only keep the inliers (set quality of outliers to zero). */
        if (covariance.array().abs().maxCoeff() < minimal_covariance) {
            for (std::size_t row = 0; row < infos->size(); ++row) {
                if (!is_inlier[row]) infos->at(row).quality = 0.0f;
            }
            return true;
        }

        /* Invert the covariance. FullPivLU is not the fastest way but
         * it gives feedback about numerical stability during inversion. */
        Eigen::FullPivLU<Eigen::Matrix3d> lu(covariance);
        if (!lu.isInvertible()) {
            return false;
        }
        covariance_inv = lu.inverse();

        /* Compute new number of inliers (all views with a gauss value above a threshold). */
        for (std::size_t row = 0; row < infos->size(); ++row) {
            Eigen::RowVector3d color = mve_to_eigen(infos->at(row).mean_color).cast<double>();
            double gauss_value = multi_gauss_unnormalized(color, var_mean, covariance_inv);
            is_inlier[row] = (gauss_value >= gauss_rejection_threshold ? 1 : 0);
        }
        /* Resize Eigen matrix accordingly and fill with new inliers. */
        inliers.resize(std::accumulate(is_inlier.begin(), is_inlier.end(), 0), Eigen::NoChange);
        for (std::size_t row = 0, inlier_row = 0; row < infos->size(); ++row) {
            if (is_inlier[row]) {
                inliers.row(inlier_row++) = mve_to_eigen(infos->at(row).mean_color).cast<double>();
            }
        }
    }

    covariance_inv *= outlier_removal_factor;
    for (ProjectedFaceInfo & info : *infos) {
        Eigen::RowVector3d color = mve_to_eigen(info.mean_color).cast<double>();
        double gauss_value = multi_gauss_unnormalized(color, var_mean, covariance_inv);
        assert(0.0 <= gauss_value && gauss_value <= 1.0);
        switch(settings.outlier_removal) {
            case NONE: return true;
            case GAUSS_DAMPING:
                info.quality *= gauss_value;
            break;
            case GAUSS_CLAMPING:
                if (gauss_value < gauss_rejection_threshold) info.quality = 0.0f;
            break;
        }
    }
    return true;
}

void
calculate_data_costs(mve::TriangleMesh::ConstPtr mesh, std::vector<TextureView> * texture_views,
    Settings const & settings, ST * data_costs) {

    mve::TriangleMesh::FaceList const & faces = mesh->get_faces();
    mve::TriangleMesh::VertexList const & vertices = mesh->get_vertices();
    mve::TriangleMesh::NormalList const & face_normals = mesh->get_face_normals();

    std::size_t const num_faces = faces.size() / 3;
    std::size_t const num_views = texture_views->size();

    CollisionModel3D* model = newCollisionModel3D(true);
    if (settings.geometric_visibility_test) {
        /* Build up acceleration structure for the visibility test. */
        ProgressCounter face_counter("\tBuilding collision model", num_faces);
        model->setTriangleNumber(num_faces);
        for (std::size_t i = 0; i < faces.size(); i += 3) {
            face_counter.progress<SIMPLE>();
            math::Vec3f v1 = vertices[faces[i]];
            math::Vec3f v2 = vertices[faces[i + 1]];
            math::Vec3f v3 = vertices[faces[i + 2]];
            model->addTriangle(*v1, *v2, *v3);
            face_counter.inc();
        }
        model->finalize();
    }
    std::vector<std::vector<ProjectedFaceInfo> > projected_face_infos(num_faces);

    ProgressCounter view_counter("\tCalculating face qualities", num_views);
    #pragma omp parallel
    {
        std::vector<std::pair<std::size_t, ProjectedFaceInfo> > projected_face_view_infos;

        #pragma omp for schedule(dynamic)
        for (std::uint16_t j = 0; j < texture_views->size(); ++j) {
            view_counter.progress<SIMPLE>();

            TextureView * texture_view = &texture_views->at(j);
            texture_view->load_image();
            texture_view->generate_validity_mask();

            if (settings.data_term == GMI) {
                texture_view->generate_gradient_magnitude();
                texture_view->erode_validity_mask();
            }

            math::Vec3f const & view_pos = texture_view->get_pos();
            math::Vec3f const & viewing_direction = texture_view->get_viewing_direction();

            for (std::size_t i = 0; i < faces.size(); i += 3) {
                std::size_t face_id = i / 3;

                math::Vec3f const & v1 = vertices[faces[i]];
                math::Vec3f const & v2 = vertices[faces[i + 1]];
                math::Vec3f const & v3 = vertices[faces[i + 2]];
                math::Vec3f const & face_normal = face_normals[face_id];
                math::Vec3f const face_center = (v1 + v2 + v3) / 3.0f;

                /* Check visibility and compute quality */

                math::Vec3f view_to_face_vec = (face_center - view_pos).normalized();
                math::Vec3f face_to_view_vec = (view_pos - face_center).normalized();

                /* Backface and basic frustum culling */
                float viewing_angle = face_to_view_vec.dot(face_normal);
                if (viewing_angle < 0.0f || viewing_direction.dot(view_to_face_vec) < 0.0f)
                    continue;

                if (std::acos(viewing_angle) > MATH_DEG2RAD(75.0f))
                    continue;

                /* Projects into the valid part of the TextureView? */
                if (!texture_view->inside(v1, v2, v3))
                    continue;

                if (settings.geometric_visibility_test) {
                    /* Viewing rays do not collide? */
                    bool visible = true;
                    math::Vec3f const * samples[] = {&v1, &v2, &v3};
                    // TODO: random monte carlo samples...

                    for (std::size_t k = 0; k < sizeof(samples) / sizeof(samples[0]); ++k) {
                        math::Vec3f vertex = *samples[k];
                        math::Vec3f dir = view_pos - vertex;
                        float const dir_length = dir.norm();
                        dir.normalize();

                        if (model->rayCollision(*vertex, *dir,  false, dir_length * 0.0001f, dir_length)) {
                            visible = false;
                            break;
                        }
                    }
                    if (!visible) continue;
                }

                ProjectedFaceInfo info = {j, 0.0f, math::Vec3f(0.0f, 0.0f, 0.0f)};

                /* Calculate quality. */
                texture_view->get_face_info(v1, v2, v3, &info, settings);

                if (info.quality == 0.0) continue;

                /* Change color space. */
                mve::image::color_rgb_to_ycbcr(*(info.mean_color));

                std::pair<std::size_t, ProjectedFaceInfo> pair(face_id, info);
                projected_face_view_infos.push_back(pair);
            }

            texture_view->release_image();
            texture_view->release_validity_mask();
            if (settings.data_term == GMI) {
                texture_view->release_gradient_magnitude();
            }
            view_counter.inc();
        }

        //std::sort(projected_face_view_infos.begin(), projected_face_view_infos.end());

        #pragma omp critical
        {
            for (std::size_t i = projected_face_view_infos.size(); 0 < i; --i) {
                std::size_t face_id = projected_face_view_infos[i - 1].first;
                ProjectedFaceInfo const & info = projected_face_view_infos[i - 1].second;
                projected_face_infos[face_id].push_back(info);
            }
            projected_face_view_infos.clear();
        }
    }

    delete model;
    model = NULL;

    ProgressCounter face_counter("\tPostprocessing face infos", num_faces);
    #pragma omp parallel for schedule(dynamic)
    for (std::size_t i = 0; i < projected_face_infos.size(); ++i) {
        face_counter.progress<SIMPLE>();

        std::vector<ProjectedFaceInfo> & infos = projected_face_infos[i];
        if (settings.outlier_removal != NONE) {
            photometric_outlier_detection(&infos, settings);

            infos.erase(std::remove_if(infos.begin(), infos.end(),
                [](ProjectedFaceInfo const & info) -> bool {return info.quality == 0.0f;}),
                infos.end());
        }
        std::sort(infos.begin(), infos.end());

        face_counter.inc();
    }

    /* Determine the function for the normlization. */
    float max_quality = 0.0f;
    for (std::size_t i = 0; i < projected_face_infos.size(); ++i)
        for (std::size_t j = 0; j < projected_face_infos[i].size(); ++j)
            max_quality = std::max(max_quality, projected_face_infos[i][j].quality);

    Histogram hist_qualities(0.0f, max_quality, 10000);
    for (std::size_t i = 0; i < projected_face_infos.size(); ++i)
        for (std::size_t j = 0; j < projected_face_infos[i].size(); ++j)
            hist_qualities.add_value(projected_face_infos[i][j].quality);

    float percentile = hist_qualities.get_approx_percentile(0.995f);

    /* Calculate the costs. */
    assert(num_faces < std::numeric_limits<std::uint32_t>::max());
    assert(num_views < std::numeric_limits<std::uint16_t>::max());
    assert(MRF_MAX_ENERGYTERM < std::numeric_limits<float>::max());
    for (std::uint32_t i = 0; i < static_cast<std::uint32_t>(projected_face_infos.size()); ++i) {
        for (std::size_t j = 0; j < projected_face_infos[i].size(); ++j) {
            ProjectedFaceInfo const & info = projected_face_infos[i][j];

            /* Clamp to percentile and normalize. */
            float normalized_quality = std::min(1.0f, info.quality / percentile);
            float data_cost = (1.0f - normalized_quality) * MRF_MAX_ENERGYTERM;
            data_costs->set_value(i, info.view_id, data_cost);
        }

        /* Ensure that all memory is freeed. */
        projected_face_infos[i] = std::vector<ProjectedFaceInfo>();
    }

    std::cout << "\tMaximum quality of a face within an image: " << max_quality << std::endl;
    std::cout << "\tClamping qualities to " << percentile << " within normalization." << std::endl;
}

TEX_NAMESPACE_END
