#include <Eigen/Core>
#include <Eigen/LU>

#include <numeric>

#include "texturing.h"
#include "SparseTable.h"
#include "Histogram.h"
#include "util.h"

bool IGNORE_LUMINANCE = false;
const int MAXVALUE = MRF_MAX_ENERGYTERM;

/** Potts model */
int
potts(int s1, int s2, int l1, int l2) {
    /* Suppress compiler warning because of unused variable. */
    (void) s1; (void) s2;
    return l1 == l2 && l1 != 0 && l2 != 0 ? 0 : 1 * MAXVALUE;
}

typedef SparseTable<std::uint32_t, std::uint16_t, int> ST;

/** Struct representing a non directed edge. */
struct Edge {
    std::size_t v1_id;
    std::size_t v2_id;

    Edge(std::size_t _v1_id, std::size_t _v2_id) {
        if(_v1_id < _v2_id) {
            v1_id = _v1_id;
            v2_id = _v2_id;
        } else {
            v1_id = _v2_id;
            v2_id = _v1_id;
        }
    }
};

/** Lexical comparison of non derected edges. */
bool
operator<(Edge const & edge1, Edge const & edge2) {
    return edge1.v1_id < edge2.v1_id || (edge1.v1_id == edge2.v1_id && edge1.v2_id < edge2.v2_id);
}

/** Setup the neighborhood of the MRF. */
void
set_neighbors(MRF * mrf, UniGraph const & graph) {
    for (std::size_t i = 0; i < graph.num_nodes(); ++i) {
        std::vector<std::size_t> adj_faces = graph.get_adj_nodes(i);
        for (std::size_t j = 0; j < adj_faces.size(); ++j) {
            std::size_t adj_face = adj_faces[j];
            /* The solver expects only one call of setNeighbours for two neighbours a and b. */
            if (i < adj_face) mrf->set_neighbors(i, adj_face);
        }
    }
}

/**
 * Dampens the quality of all views in which the face's projection
 * has a much different color than in the majority of views.
 *
 * @param infos Contains information about one face seen from several views
 * @param conf The program configuration.
 */
void
photometric_outlier_detection(std::vector<ProjectedFaceInfo> & infos, Arguments const & conf) {
    if (conf.outlier_removal == NONE) return;
    if (infos.size() == 0) return;

    /* Configuration variables. */

    double const gauss_rejection_threshold = 6e-3;
    /* If all covariances drop below this we stop outlier detection. */
    double const minimal_covariance = 5e-4;

    int const outlier_detection_iterations = 10;
    int const minimal_num_inliers = 4;

    float outlier_removal_factor = std::numeric_limits<float>::signaling_NaN();
    switch (conf.outlier_removal) {
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
        switch(conf.outlier_removal) {
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

/** Calculate data costs for the given data term (conf.data_term). */
ST
calculate_data_costs(mve::TriangleMesh::ConstPtr mesh, std::vector<TextureView> const & texture_views,
    Arguments const & conf) {

    mve::TriangleMesh::FaceList const & faces = mesh->get_faces();
    mve::TriangleMesh::VertexList const & vertices = mesh->get_vertices();
    mve::TriangleMesh::NormalList const & face_normals = mesh->get_face_normals();

    std::size_t const num_faces = faces.size() / 3;
    std::size_t const num_views = texture_views.size();

    CollisionModel3D* model = newCollisionModel3D(true);
    if (conf.geometric_visibility_test) {
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

            /* Check visibility and compute quality of each face in each texture view. */
            for (std::uint16_t j = 0; j < texture_views.size(); ++j) {
                math::Vec3f const & v1 = vertices[faces[i]];
                math::Vec3f const & v2 = vertices[faces[i + 1]];
                math::Vec3f const & v3 = vertices[faces[i + 2]];

                math::Vec3f const & view_pos = texture_views[j].get_pos();

                math::Vec3f const & face_normal = face_normals[face_id];
                math::Vec3f face_center = (v1 + v2 + v3) / 3;
                math::Vec3f view_to_face_vec = (face_center - view_pos).normalized();
                math::Vec3f face_to_view_vec = (view_pos - face_center).normalized();
                math::Vec3f const & viewing_direction = texture_views[j].get_viewing_direction();

                /* Backface culling */
                float viewing_angle = face_to_view_vec.dot(face_normal);
                if (viewing_angle < 0.0f || viewing_direction.dot(view_to_face_vec) < 0.0f)
                    continue;

                /* Projects into the valid part of the TextureView? */
                if (!texture_views[j].inside(v1, v2, v3))
                    continue;

                if (conf.geometric_visibility_test) {
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
                    if (!visible)
                        continue;
                }

                ProjectedFaceInfo info = {j, 0.0f, math::Vec3f(0.0f, 0.0f, 0.0f)};

                /* Calculate quality. */
                switch (conf.data_term) {
                    case AREA:
                        texture_views[j].get_face_info<AREA>(v1, v2, v3, &info);
                    break;
                    case GMI:
                        texture_views[j].get_face_info<GMI>(v1, v2, v3, &info);
                    break;
                }

                if (info.quality == 0.0f)
                    continue;

                /* Change color space. */
                mve::image::color_rgb_to_ycbcr(*(info.mean_color));

                infos.push_back(info);
            }

            photometric_outlier_detection(infos, conf);

            reduced_projected_face_infos[face_id].reserve(infos.size());
            for (ProjectedFaceInfo const & info : infos)
                reduced_projected_face_infos[face_id].push_back(reduce(info));

            face_counter.inc();
        }
    }

	delete model;
	model = NULL;

    /* Determine the function for the quantization. */
    float max_quality = 0.0f;
    for (std::size_t i = 0; i < reduced_projected_face_infos.size(); ++i)
        for (std::size_t j = 0; j < reduced_projected_face_infos[i].size(); ++j)
            max_quality = std::max(max_quality, reduced_projected_face_infos[i][j].quality);

    Histogram hist_qualities(0.0f, max_quality, 1000);
    for (std::size_t i = 0; i < reduced_projected_face_infos.size(); ++i)
        for (std::size_t j = 0; j < reduced_projected_face_infos[i].size(); ++j)
            hist_qualities.add_value(reduced_projected_face_infos[i][j].quality);

    if (conf.write_data_term_histograms)
        hist_qualities.save_to_file(conf.out_prefix + "_hist_qualities.csv");

    float permille_999 = hist_qualities.get_approximate_permille(0.999f);

    Histogram hist_nqualities(0.0f, 1.0f, 1000);

    /* Calculate the costs. */
    assert(num_faces < std::numeric_limits<std::uint32_t>::max());
    assert(num_views < std::numeric_limits<std::uint16_t>::max());
    assert(MAXVALUE < std::numeric_limits<int>::max());
    ST data_costs(num_faces, num_views);
    for (std::uint32_t i = 0; i < static_cast<std::uint32_t>(reduced_projected_face_infos.size()); ++i) {
        while(!reduced_projected_face_infos[i].empty()) {
            ReducedProjectedFaceInfo info = reduced_projected_face_infos[i].back();
            reduced_projected_face_infos[i].pop_back();

            /* Clamp to percentile and normalize. */
            float normalized_quality = std::min(1.0f, info.quality / permille_999);
            hist_nqualities.add_value(normalized_quality);
            int data_cost = (1.0f - normalized_quality) * MAXVALUE;
            data_costs.set_value(i, info.view_id, data_cost);
        }

        /* Ensure that all memory is freeed. */
        reduced_projected_face_infos[i].clear();
        reduced_projected_face_infos[i].shrink_to_fit();
    }
    if (conf.write_data_term_histograms)
        hist_nqualities.save_to_file(conf.out_prefix + "_hist_nqualities.csv");

    std::cout << "\tMaximum quality of a face within an image: " << max_quality << std::endl;
    std::cout << "\tClamping qualities to " << permille_999 << " within discretization." << std::endl;
    return data_costs;
}


/** Set the data costs of the MRF. */
void
set_data_costs(MRF * mrf, ST const & data_costs){
    /* Set data costs for all labels except label 0 (undefined) */
    for (std::size_t i = 0; i < data_costs.rows(); i++) {
        ST::Row const & data_costs_for_label = data_costs.row(i);

        std::vector<MRF::SparseDataCost> costs(data_costs_for_label.size());
        for(std::size_t j = 0; j < costs.size(); j++) {
            costs[j].site = data_costs_for_label[j].first;
            costs[j].cost = data_costs_for_label[j].second;
        }


        int label = i + 1;
        mrf->set_data_costs(label, &costs);
    }

    /* Set costs for undefined label */
    std::vector<MRF::SparseDataCost> costs(data_costs.cols());
    for (std::size_t i = 0; i < costs.size(); i++) {
        costs[i].site = i;
        costs[i].cost = MRF_MAX_ENERGYTERM;
    }
    mrf->set_data_costs(0, &costs);
}

/** Remove all edges of nodes which corresponding face has not been seen in any texture view. */
void
isolate_unseen_faces(UniGraph * graph, ST const & data_costs) {
    int num_unseen_faces = 0;
    for (std::uint32_t i = 0; i < data_costs.cols(); i++) {
        ST::Column const & data_costs_for_face = data_costs.col(i);

        if (data_costs_for_face.size() == 0) {
            num_unseen_faces++;

            std::vector<std::size_t> const & adj_nodes = graph->get_adj_nodes(i);
            for (std::size_t j = 0; j < adj_nodes.size(); j++)
                graph->remove_edge(i, adj_nodes[j]);
        }

    }
    std::cout << "\t" << num_unseen_faces << " faces have not been seen by a view." << std::endl;
}

void
build_mrf(MRF * mrf, mve::TriangleMesh::ConstPtr mesh, std::vector<TextureView>
    const & texture_views, UniGraph const & graph, Arguments const & conf) {
    std::size_t const num_faces = mesh->get_faces().size() / 3;

    ST data_costs;
    if (conf.data_cost_file.empty()) {
        data_costs = calculate_data_costs(mesh, texture_views, conf);

        if (conf.write_intermediate_results) {
            std::cout << "\tWriting data cost file... " << std::flush;
                ST::save_to_file(data_costs, conf.out_prefix + "_data_costs.spt");
            std::cout << "done." << std::endl;
        }
    } else {
        std::cout << "\tLoading data cost file... " << std::flush;
        data_costs = ST::load_from_file(conf.data_cost_file);
        if (data_costs.rows() != texture_views.size() || data_costs.cols() != num_faces) {
            std::cout << "failed!" << std::endl;
            std::cerr << "Wrong datacost file for this mesh/scene combination... aborting!" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        std::cout << "done." << std::endl;
    }

    {
        UniGraph mod_graph(graph); //TODO temporary waste of memory...
        isolate_unseen_faces(&mod_graph, data_costs);

        /* Set neighbors must be called prior to set_data_costs (LBP). */
        set_neighbors(mrf, mod_graph);
    }

    set_data_costs(mrf, data_costs);

    switch (conf.smoothness_term) {
        case EDI:
            std::cerr << "\tThe smoothness term EDI has been removed - falling back to POTTS" << std::endl;
        case POTTS:
            mrf->set_smooth_cost(*potts);
        break;
    }
}
