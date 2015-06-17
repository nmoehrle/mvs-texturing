#include <Eigen/Core>
#include <Eigen/LU>

#include <numeric>

#include "texturing.h"
#include "SparseTable.h"
#include "Histogram.h"

TEX_NAMESPACE_BEGIN

/**
 * Dampens the quality of all views in which the face's projection
 * has a much different color than in the majority of views.
 *
 * @param infos Contains information about one face seen from several views
 * @param conf The program configuration.
 */
void
photometric_outlier_detection(std::vector<ProjectedFaceInfo> & infos, Settings const & settings) {
    if (settings.outlier_removal == NONE) return;
    if (infos.size() == 0) return;

    /* Configuration variables. */

    double const gauss_rejection_threshold = 6e-3;
    /* If all covariances drop below this we stop outlier detection. */
    double const minimal_covariance = 5e-4;

    int const outlier_detection_iterations = 10;
    int const minimal_num_inliers = 4;

    float outlier_removal_factor = std::numeric_limits<float>::signaling_NaN();
    switch (settings.outlier_removal) {
        case NONE: return;
        case GAUSS_CLAMPING:
            outlier_removal_factor = 1.0f;
        break;
        case GAUSS_DAMPING:
            outlier_removal_factor = 0.2f;
        break;
    }

    Eigen::MatrixX3d inliers(infos.size(), 3);
    std::vector<std::uint32_t> is_inlier(infos.size());
    for (std::size_t row = 0; row < infos.size(); ++row) {
        inliers.row(row) = mve_to_eigen(infos.at(row).mean_color).cast<double>();
        is_inlier.at(row) = 1;
    }

    Eigen::RowVector3d var_mean;
    Eigen::Matrix3d covariance;
    Eigen::Matrix3d covariance_inv;

    /* This variable indicates whether something went wrong in the outlier detection
     * (number of inliers below threshold or inverting the covariance goes wrong).
     * In this case no outlier removal is done afterwards. */
    bool bad_outlier_detection = false;

    for (int outlier_detection_it = 0; outlier_detection_it < outlier_detection_iterations; ++outlier_detection_it) {

        if (inliers.rows() < minimal_num_inliers) {
            bad_outlier_detection = true;
            break;
        }

        /* Calculate the inliers' mean color and color covariance. */
        var_mean = inliers.colwise().mean();
        Eigen::MatrixX3d centered = inliers.rowwise() - var_mean;
        covariance = (centered.adjoint() * centered) / double(inliers.rows() - 1);

        /* If all covariances are very small we stop outlier detection
         * (this does not mean it went wrong so we don't set the variable). */
        if (covariance.array().abs().maxCoeff() < minimal_covariance)
            break;

        /* Invert the covariance. FullPivLU is not the fastest way but
         * it gives feedback about numerical stability during inversion. */
        Eigen::FullPivLU<Eigen::Matrix3d> lu(covariance);
        if (!lu.isInvertible()) {
            bad_outlier_detection = true;
            break;
        }
        covariance_inv = lu.inverse();

        /* Compute new number of inliers (all views with a gauss value above a threshold). */
        for (std::size_t row = 0; row < infos.size(); ++row) {
            Eigen::RowVector3d color = mve_to_eigen(infos.at(row).mean_color).cast<double>();
            double gauss_value = multi_gauss_unnormalized(color, var_mean, covariance_inv);
            is_inlier.at(row) = (gauss_value >= gauss_rejection_threshold ? 1 : 0);
        }
        /* Resize Eigen matrix accordingly and fill with new inliers. */
        inliers.resize(std::accumulate(is_inlier.begin(), is_inlier.end(), 0), Eigen::NoChange);
        for (std::size_t row = 0, inlier_row = 0; row < infos.size(); ++row) {
            if (is_inlier.at(row)) {
                inliers.row(inlier_row) = mve_to_eigen(infos.at(row).mean_color).cast<double>();
                inlier_row++;
            }
        }
    }

    /* If something went wrong during outlier detection we don't do any removal. */
    if(bad_outlier_detection)
        return;

    /* If the covariance is very very small we only keep the inliers
     * no matter which detection method has been chosen. */
    if (covariance.array().abs().maxCoeff() < minimal_covariance) {
        for (std::size_t row = 0; row < infos.size(); ++row)
            if (is_inlier.at(row)) infos.at(row).quality = 0.0f;
        return;
    }

    covariance_inv *= outlier_removal_factor;
    for (ProjectedFaceInfo & info : infos) {
        Eigen::RowVector3d color = mve_to_eigen(info.mean_color).cast<double>();
        double gauss_value = multi_gauss_unnormalized(color, var_mean, covariance_inv);
        assert(0.0 <= gauss_value && gauss_value <= 1.0);
        switch(settings.outlier_removal) {
            case NONE: return;
            case GAUSS_DAMPING:
                info.quality *= gauss_value;
            break;
            case GAUSS_CLAMPING:
                if (gauss_value < gauss_rejection_threshold) info.quality = 0;
            break;
        }
    }
}

ST
calculate_data_costs(mve::TriangleMesh::ConstPtr mesh, std::vector<TextureView> & texture_views,
    Settings const & settings) {

    mve::TriangleMesh::FaceList const & faces = mesh->get_faces();
    mve::TriangleMesh::VertexList const & vertices = mesh->get_vertices();
    mve::TriangleMesh::NormalList const & face_normals = mesh->get_face_normals();

    std::size_t const num_faces = faces.size() / 3;
    std::size_t const num_views = texture_views.size();

    ProgressCounter view_counter("\tGenerating validity masks", num_views);
    #pragma omp parallel for
    for (std::size_t i = 0; i < num_views; ++i) {
        view_counter.progress<SIMPLE>();
        texture_views[i].generate_validity_mask();
        view_counter.inc();
    }

    if (settings.data_term == GMI) {
        view_counter.reset("\tCalculating gradient magnitude");
        #pragma omp parallel for
        for (std::size_t i = 0; i < num_views; ++i) {
            view_counter.progress<SIMPLE>();
            texture_views[i].generate_gradient_magnitude();
            texture_views[i].erode_validity_mask();
            view_counter.inc();
        }
    }

    CollisionModel3D* model = newCollisionModel3D(true);
    if (settings.geometric_visibility_test) {
        /* Build up acceleration structure for the visibility test. */
        ProgressCounter face_counter("\tBuilding collision model for visibility tests", num_faces);
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
    std::vector<std::vector<ReducedProjectedFaceInfo> > reduced_projected_face_infos(num_faces);

    ProgressCounter face_counter("\tCalculating face qualities", num_faces);
    #pragma omp parallel
    {
        std::vector<ProjectedFaceInfo> infos;
        infos.reserve(num_views);

        #pragma omp for schedule(dynamic)
        for (std::size_t i = 0; i < faces.size(); i += 3) {
            std::size_t face_id = i / 3;

            face_counter.progress<ETA>();
            infos.clear();
            infos.reserve(num_views);

            math::Vec3f const & v1 = vertices[faces[i]];
            math::Vec3f const & v2 = vertices[faces[i + 1]];
            math::Vec3f const & v3 = vertices[faces[i + 2]];
            math::Vec3f const & face_normal = face_normals[face_id];
            math::Vec3f const face_center = (v1 + v2 + v3) / 3;

            /* Check visibility and compute quality of each face in each texture view. */
            for (std::uint16_t j = 0; j < texture_views.size(); ++j) {
                math::Vec3f const & view_pos = texture_views[j].get_pos();
                math::Vec3f const & viewing_direction = texture_views[j].get_viewing_direction();

                math::Vec3f view_to_face_vec = (face_center - view_pos).normalized();
                math::Vec3f face_to_view_vec = (view_pos - face_center).normalized();

                /* Backface culling */
                float viewing_angle = face_to_view_vec.dot(face_normal);
                if (viewing_angle < 0.0f || viewing_direction.dot(view_to_face_vec) < 0.0f)
                    continue;

                /* Projects into the valid part of the TextureView? */
                if (!texture_views[j].inside(v1, v2, v3))
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
                texture_views[j].get_face_info(v1, v2, v3, &info, settings);

                if (info.quality == 0.0) continue;

                /* Change color space. */
                mve::image::color_rgb_to_ycbcr(*(info.mean_color));

                infos.push_back(info);
            }

            photometric_outlier_detection(infos, settings);

            reduced_projected_face_infos[face_id].reserve(infos.size());
            for (ProjectedFaceInfo const & info : infos) {
                reduced_projected_face_infos[face_id].push_back(reduce(info));
            }

            face_counter.inc();
        }
    }

    delete model;
    model = NULL;

    /* Determine the function for the normlization. */
    float max_quality = 0.0f;
    for (std::size_t i = 0; i < reduced_projected_face_infos.size(); ++i)
        for (std::size_t j = 0; j < reduced_projected_face_infos[i].size(); ++j)
            max_quality = std::max(max_quality, reduced_projected_face_infos[i][j].quality);

    Histogram hist_qualities(0.0f, max_quality, 1000);
    for (std::size_t i = 0; i < reduced_projected_face_infos.size(); ++i)
        for (std::size_t j = 0; j < reduced_projected_face_infos[i].size(); ++j)
            hist_qualities.add_value(reduced_projected_face_infos[i][j].quality);

    float permille_999 = hist_qualities.get_approximate_permille(0.999f);

    /* Calculate the costs. */
    assert(num_faces < std::numeric_limits<std::uint32_t>::max());
    assert(num_views < std::numeric_limits<std::uint16_t>::max());
    assert(MRF_MAX_ENERGYTERM < std::numeric_limits<float>::max());
    ST data_costs(num_faces, num_views);
    for (std::uint32_t i = 0; i < static_cast<std::uint32_t>(reduced_projected_face_infos.size()); ++i) {
        while(!reduced_projected_face_infos[i].empty()) {
            ReducedProjectedFaceInfo info = reduced_projected_face_infos[i].back();
            reduced_projected_face_infos[i].pop_back();

            /* Clamp to percentile and normalize. */
            float normalized_quality = std::min(1.0f, info.quality / permille_999);
            float data_cost = (1.0f - normalized_quality) * MRF_MAX_ENERGYTERM;
            data_costs.set_value(i, info.view_id, data_cost);
        }

        /* Ensure that all memory is freeed. */
        reduced_projected_face_infos[i].clear();
        reduced_projected_face_infos[i].shrink_to_fit();
    }

    std::cout << "\tMaximum quality of a face within an image: " << max_quality << std::endl;
    std::cout << "\tClamping qualities to " << permille_999 << " within normalization." << std::endl;

    /* Release superfluous embeddings. */
    for (TextureView & texture_view : texture_views) {
        texture_view.release_validity_mask();
        if (settings.data_term == GMI) {
            texture_view.release_gradient_magnitude();
        }
    }

    return data_costs;
}

TEX_NAMESPACE_END