/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <map>
#include <set>

#include <util/timer.h>
#include <math/accum.h>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

#include "texturing.h"
#include "seam_leveling.h"
#include "progress_counter.h"

TEX_NAMESPACE_BEGIN

typedef Eigen::SparseMatrix<float> SpMat;

math::Vec3f
sample_edge(TexturePatch::ConstPtr texture_patch, math::Vec2f p1, math::Vec2f p2) {
    math::Vec2f p12 = p2 - p1;

    std::size_t num_samples = std::max(p12.norm(), 1.0f) * 2.0f;

    math::Accum<math::Vec3f> color_accum(math::Vec3f(0.0f));
    /* Sample the edge with linear weights. */
    for (std::size_t s = 0; s < num_samples; ++s) {
        float fraction = static_cast<float>(s) / (num_samples - 1);
        math::Vec2f sample_point = p1 + p12 * fraction;
        math::Vec3f color(texture_patch->get_pixel_value(sample_point));

        color_accum.add(color, 1.0f - fraction);
    }

    return color_accum.normalized();
}


void
find_seam_edges_for_vertex_label_combination(UniGraph const & graph, mve::TriangleMesh::ConstPtr & mesh,
    mve::MeshInfo const & mesh_info, std::size_t vertex, std::size_t label1, std::size_t label2,
    std::vector<MeshEdge> * seam_edges) {

    assert(label1 != 0 && label2 != 0 && label1 < label2);

    mve::TriangleMesh::VertexList const & vertices = mesh->get_vertices();

    std::vector<std::size_t> const & adj_verts = mesh_info[vertex].verts;
    for (std::size_t i = 0; i < adj_verts.size(); ++i) {
        std::size_t adj_vertex = adj_verts[i];
        if (vertex == adj_vertex) continue;

        std::vector<std::size_t> edge_faces;
        mesh_info.get_faces_for_edge(vertex, adj_vertex, &edge_faces);

        for (std::size_t j = 0; j < edge_faces.size(); ++j) {
            for(std::size_t k = j + 1; k < edge_faces.size(); ++k) {

                std::size_t face_label1 = graph.get_label(edge_faces[j]);
                std::size_t face_label2 = graph.get_label(edge_faces[k]);
                if (!(face_label1 < face_label2)) std::swap(face_label1, face_label2);

                if (face_label1 != label1 || face_label2 != label2) continue;

                math::Vec3f v1 = vertices[vertex];
                math::Vec3f v2 = vertices[adj_vertex];
                float length = (v2 - v1).norm();

                /* Ignore zero length edges. */
                if (length == 0.0f) continue;

                MeshEdge seam_edge = {vertex, adj_vertex};
                seam_edges->push_back(seam_edge);
            }
        }
    }
}

math::Vec3f
calculate_difference(VertexProjectionInfos const & vertex_projection_infos,
    mve::TriangleMesh::ConstPtr & mesh, std::vector<TexturePatch::Ptr> const &  texture_patches,
    std::vector<MeshEdge> const & seam_edges, int label1, int label2) {

    assert(label1 != 0 && label2 != 0 && label1 < label2);
    assert(!seam_edges.empty());

    mve::TriangleMesh::VertexList const & vertices = mesh->get_vertices();

    math::Accum<math::Vec3f> color1_accum(math::Vec3f(0.0f));
    math::Accum<math::Vec3f> color2_accum(math::Vec3f(0.0f));

    for (MeshEdge const & seam_edge : seam_edges) {
        math::Vec3f v1 = vertices[seam_edge.v1];
        math::Vec3f v2 = vertices[seam_edge.v2];
        float length = (v2 - v1).norm();

        assert(length != 0.0f);

        std::vector<EdgeProjectionInfo> edge_projection_infos;
        find_mesh_edge_projections(vertex_projection_infos, seam_edge, &edge_projection_infos);

        std::size_t num_samples = 0;

        for (EdgeProjectionInfo const & edge_projection_info : edge_projection_infos) {
            TexturePatch::Ptr texture_patch = texture_patches[edge_projection_info.texture_patch_id];
            const int texture_patch_label = texture_patch->get_label();
            if (texture_patch_label == label1 || texture_patch_label == label2) {
                if (texture_patch_label == label1)
                    color1_accum.add(sample_edge(texture_patch, edge_projection_info.p1, edge_projection_info.p2), length);

                if (texture_patch_label == label2)
                    color2_accum.add(sample_edge(texture_patch, edge_projection_info.p1, edge_projection_info.p2), length);

                num_samples++;
            }
        }
        assert(num_samples == 2);
    }

    math::Vec3f color1 = color1_accum.normalized();
    math::Vec3f color2 = color2_accum.normalized();

    /* The order is essential. */
    math::Vec3f difference = color2 - color1;

    assert(!std::isnan(difference[0]));
    assert(!std::isnan(difference[1]));
    assert(!std::isnan(difference[2]));

    return difference;
}

void
global_seam_leveling(UniGraph const & graph, mve::TriangleMesh::ConstPtr mesh,
    mve::MeshInfo const & mesh_info,
    std::vector<std::vector<VertexProjectionInfo> > const & vertex_projection_infos,
    std::vector<TexturePatch::Ptr> * texture_patches) {

    mve::TriangleMesh::VertexList const & vertices = mesh->get_vertices();
    std::size_t const num_vertices = vertices.size();

    std::cout << "\tCreate matrices for optimization... " << std::flush;
    std::vector<std::map<std::size_t, std::size_t> > vertlabel2row;
    vertlabel2row.resize(num_vertices);

    std::vector<std::vector<std::size_t> > labels;
    labels.resize(num_vertices);

    /* Assign each vertex for each label a new index(row) within the solution vector x. */
    std::size_t x_row = 0;
    for (std::size_t i = 0; i < num_vertices; ++i) {
        std::set<std::size_t> label_set;

        std::vector<std::size_t> faces = mesh_info[i].faces;
        std::set<std::size_t>::iterator it = label_set.begin();
        for (std::size_t j = 0; j < faces.size(); ++j) {
            std::size_t label = graph.get_label(faces[j]);
            label_set.insert(it, label);
        }

        for (it = label_set.begin(); it != label_set.end(); ++it) {
            std::size_t label = *it;
            if (label == 0) continue;
            vertlabel2row[i][label] = x_row;
            labels[i].push_back(label);
            ++x_row;
        }
    }
    std::size_t x_rows = x_row;
    assert(x_rows < static_cast<std::size_t>(std::numeric_limits<int>::max()));

    float const lambda = 0.1f;

    /* Fill the Tikhonov matrix Gamma(regularization constraints). */
    std::size_t Gamma_row = 0;
    std::vector<Eigen::Triplet<float, int> > coefficients_Gamma;
    coefficients_Gamma.reserve(2 * num_vertices);
    for (std::size_t i = 0; i < num_vertices; ++i) {
        for (std::size_t j = 0; j < labels[i].size(); ++j) {
            std::vector<std::size_t> const & adj_verts = mesh_info[i].verts;
            for (std::size_t k = 0; k < adj_verts.size(); ++k) {
                std::size_t adj_vertex = adj_verts[k];
                for (std::size_t l = 0; l < labels[adj_vertex].size(); ++l) {
                    std::size_t label = labels[i][j];
                    std::size_t adj_vertex_label = labels[adj_vertex][l];
                    if (i < adj_vertex && label == adj_vertex_label) {
                        Eigen::Triplet<float, int> t1(Gamma_row, vertlabel2row[i][label], lambda);
                        Eigen::Triplet<float, int> t2(Gamma_row, vertlabel2row[adj_vertex][adj_vertex_label], -lambda);
                        coefficients_Gamma.push_back(t1);
                        coefficients_Gamma.push_back(t2);
                        Gamma_row++;
                    }
                }
            }
        }
    }
    std::size_t Gamma_rows = Gamma_row;
    assert(Gamma_rows < static_cast<std::size_t>(std::numeric_limits<int>::max()));

    SpMat Gamma(Gamma_rows, x_rows);
    Gamma.setFromTriplets(coefficients_Gamma.begin(), coefficients_Gamma.end());

    /* Fill the matrix A and the coefficients for the Vector b of the linear equation system. */
    std::vector<Eigen::Triplet<float, int> > coefficients_A;
    std::vector<math::Vec3f> coefficients_b;
    std::size_t A_row = 0;
    for (std::size_t i = 0; i < num_vertices; ++i) {
        for (std::size_t j = 0; j < labels[i].size(); ++j) {
            for (std::size_t k = 0; k < labels[i].size(); ++k) {
                std::size_t label1 = labels[i][j];
                std::size_t label2 = labels[i][k];
                if (label1 < label2) {

                    std::vector<MeshEdge> seam_edges;
                    find_seam_edges_for_vertex_label_combination(graph, mesh, mesh_info, i, label1, label2, &seam_edges);

                    if (seam_edges.empty()) continue;

                    Eigen::Triplet<float, int> t1(A_row, vertlabel2row[i][label1], 1.0f);
                    Eigen::Triplet<float, int> t2(A_row, vertlabel2row[i][label2], -1.0f);
                    coefficients_A.push_back(t1);
                    coefficients_A.push_back(t2);

                    coefficients_b.push_back(calculate_difference(vertex_projection_infos, mesh, *texture_patches, seam_edges, label1, label2));

                    ++A_row;
                }
            }
        }
    }

    std::size_t A_rows = A_row;
    assert(A_rows < static_cast<std::size_t>(std::numeric_limits<int>::max()));

    SpMat A(A_rows, x_rows);
    A.setFromTriplets(coefficients_A.begin(), coefficients_A.end());

    SpMat I(x_rows, x_rows);
    I.setIdentity();

    SpMat Lhs = A.transpose() * A + Gamma.transpose() * Gamma + I * 0.0001f;

    /* Only keep lower triangle (CG only uses the lower),
     * prune the rest and compress matrix. */
    Lhs.prune([](const int& row, const int& col, const float& value) -> bool {
            return col <= row && value != 0.0f;
        }); // value != 0.0f is only to suppress a compiler warning

    std::vector<std::map<std::size_t, math::Vec3f> > adjust_values(num_vertices);
    std::cout << " done." << std::endl;
    std::cout << "\tLhs dimensionality: " << Lhs.rows() << " x " << Lhs.cols() << std::endl;

    util::WallTimer timer;
    std::cout << "\tCalculating adjustments:"<< std::endl;
    #pragma omp parallel for
    for (std::size_t channel = 0; channel < 3; ++channel) {
        /* Prepare solver. */
        Eigen::ConjugateGradient<SpMat, Eigen::Lower> cg;
        cg.setMaxIterations(1000);
        cg.setTolerance(0.0001);
        cg.compute(Lhs);

        /* Prepare right hand side. */
        Eigen::VectorXf b(A_rows);
        for (std::size_t i = 0; i < coefficients_b.size(); ++i) {
            b[i] = coefficients_b[i][channel];
        }
        Eigen::VectorXf Rhs = SpMat(A.transpose()) * b;

        /* Solve for x. */
        Eigen::VectorXf x(x_rows);
        x = cg.solve(Rhs);

        /* Subtract mean because system is underconstrained and we seek the solution with minimal adjustments. */
        x = x.array() - x.mean();

        #pragma omp critical
        std::cout << "\t\tColor channel " << channel << ": CG took "
            << cg.iterations() << " iterations. Residual is " << cg.error() << std::endl;

        #pragma omp critical
        for (std::size_t i = 0; i < num_vertices; ++i) {
            for (std::size_t j = 0; j < labels[i].size(); ++j) {
                std::size_t label = labels[i][j];
                adjust_values[i][label][channel] = x[vertlabel2row[i][label]];
            }
        }
    }
    std::cout << "\t\tTook " << timer.get_elapsed_sec() << " seconds" << std::endl;

    mve::TriangleMesh::FaceList const & mesh_faces = mesh->get_faces();

    ProgressCounter texture_patch_counter("\tAdjusting texture patches", texture_patches->size());
    #pragma omp parallel for schedule(dynamic)
    for (std::size_t i = 0; i < texture_patches->size(); ++i) {
        texture_patch_counter.progress<SIMPLE>();

        TexturePatch::Ptr texture_patch = texture_patches->at(i);

        int label = texture_patch->get_label();
        std::vector<std::size_t> const & faces = texture_patch->get_faces();
        std::vector<math::Vec3f> patch_adjust_values(faces.size() * 3, math::Vec3f(0.0f));

        /* Only adjust texture_patches originating form input images. */
        if (label == 0) {
            texture_patch->adjust_colors(patch_adjust_values);
            texture_patch_counter.inc();
            continue;
        };

        for (std::size_t j = 0; j < faces.size(); ++j) {
            for (std::size_t k = 0; k < 3; ++k) {
                std::size_t face_pos = faces[j] * 3 + k;
                std::size_t vertex = mesh_faces[face_pos];
                patch_adjust_values[j * 3 + k] = adjust_values[vertex].find(label)->second;
            }
        }

        texture_patch->adjust_colors(patch_adjust_values);
        texture_patch_counter.inc();
    }
}

TEX_NAMESPACE_END
