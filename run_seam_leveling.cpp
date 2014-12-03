#include "texturing.h"
#include <Eigen/Sparse>
#include <map>
#include <set>
#include "util/timer.h"

#define STRIP_SIZE 20

typedef Eigen::SparseMatrix<float> SpMat;

struct ProjectedEdgeInfo {
    std::size_t texture_patch_id;
    math::Vec2f p1;
    math::Vec2f p2;

    bool operator<(ProjectedEdgeInfo const & projected_edge_info) const {
        return texture_patch_id < projected_edge_info.texture_patch_id;
    }
};

struct MeshEdge {
    std::size_t v1;
    std::size_t v2;
};

void
find_seam_edges(UniGraph const & graph, mve::TriangleMesh::ConstPtr mesh,
    std::vector<MeshEdge> * seam_edges) {
    mve::TriangleMesh::FaceList const & faces = mesh->get_faces();

    seam_edges->clear();

    // Is it possible that a single edge is part of more than three faces whichs' label is non zero???

    for (std::size_t node = 0; node < graph.num_nodes(); ++node) {
        std::vector<std::size_t> const & adj_nodes = graph.get_adj_nodes(node);
        for (std::size_t adj_node : adj_nodes) {
            /* Add each edge only once. */
            if (node > adj_node) continue;

            int label1 = graph.get_label(node);
            int label2 = graph.get_label(adj_node);
            /* Add only seam edges. */
            if (label1 == 0 || label2 == 0 || label1 == label2) continue;

            /* Find shared edge of the faces. */
            std::vector<std::size_t> shared_edge;
            for (int i = 0; i < 3; ++i){
                std::size_t v1 = faces[3 * node + i];

                for (int j = 0; j < 3; j++){
                    std::size_t v2 = faces[3 * adj_node + j];

                    if (v1 == v2) shared_edge.push_back(v1);
                }
            }

            assert(shared_edge.size() == 2);
            std::size_t v1 = shared_edge[0];
            std::size_t v2 = shared_edge[1];

            assert(v1 != v2);
            if (v1 > v2) std::swap(v1, v2);

            MeshEdge seam_edge = {v1, v2};
            seam_edges->push_back(seam_edge);
        }
    }
}

void
find_mesh_edge_projections(std::vector<std::vector<VertexProjectionInfo> > const & vertex_projection_infos,
    MeshEdge mesh_edge, std::vector<ProjectedEdgeInfo> * projected_edge_infos) {
    std::vector<VertexProjectionInfo> const & v1_projection_infos = vertex_projection_infos[mesh_edge.v1];
    std::vector<VertexProjectionInfo> const & v2_projection_infos = vertex_projection_infos[mesh_edge.v2];

    /* Use a set to eliminate duplicates which may occur if the mesh is degenerated. */
    std::set<ProjectedEdgeInfo> projected_edge_infos_set;

    for (VertexProjectionInfo v1_projection_info : v1_projection_infos) {
        for (VertexProjectionInfo v2_projection_info : v2_projection_infos) {
            if (v1_projection_info.texture_patch_id != v2_projection_info.texture_patch_id) continue;

            for (std::size_t face_id1 : v1_projection_info.faces) {
                for (std::size_t face_id2 : v2_projection_info.faces) {
                    if(face_id1 != face_id2) continue;

                    std::size_t texture_patch_id = v1_projection_info.texture_patch_id;
                    math::Vec2f p1 = v1_projection_info.projection;
                    math::Vec2f p2 = v2_projection_info.projection;

                    ProjectedEdgeInfo projected_edge_info = {texture_patch_id, p1, p2};
                    projected_edge_infos_set.insert(projected_edge_info);
                }
            }
        }
    }

    projected_edge_infos->insert(projected_edge_infos->end(), projected_edge_infos_set.begin(), projected_edge_infos_set.end());
}

math::Vec3f
sample_edge(TexturePatch const & texture_patch, math::Vec2f p1, math::Vec2f p2) {
    math::Vec2f p12 = p2 - p1;

    std::size_t num_samples = std::max(p12.norm(), 1.0f) * 2.0f;

    math::Accum<math::Vec3f> color_accum(math::Vec3f(0.0f));
    /* Sample the edge with linear weights. */
    for (std::size_t s = 0; s < num_samples; ++s) {
        float fraction = static_cast<float>(s) / (num_samples - 1);
        math::Vec2f sample_point = p1 + p12 * fraction;
        math::Vec3f color(texture_patch.get_pixel_value(sample_point));

        color_accum.add(color, 1.0f - fraction);
    }

    return color_accum.normalized();
}


void
find_seam_edges_for_vertex_label_combination(UniGraph const & graph, mve::TriangleMesh::ConstPtr & mesh,
    mve::VertexInfoList::ConstPtr vertex_infos, std::size_t vertex, std::size_t label1, std::size_t label2,
    std::vector<MeshEdge> * seam_edges) {

    assert(label1 != 0 && label2 != 0 && label1 < label2);

    mve::TriangleMesh::VertexList const & vertices = mesh->get_vertices();

    std::vector<std::size_t> const & adj_verts = vertex_infos->at(vertex).verts;
    for (std::size_t i = 0; i < adj_verts.size(); ++i) {
        std::size_t adj_vertex = adj_verts[i];
        if (vertex == adj_vertex) continue;

        std::vector<std::size_t> edge_faces;
        vertex_infos->get_faces_for_edge(vertex, adj_vertex, &edge_faces);

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
calculate_difference(std::vector<std::vector<VertexProjectionInfo> > const & vertex_projection_infos,
    mve::TriangleMesh::ConstPtr & mesh, std::vector<TexturePatch> const &  texture_patches,
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

        std::vector<ProjectedEdgeInfo> projected_edge_infos;
        find_mesh_edge_projections(vertex_projection_infos, seam_edge, &projected_edge_infos);

        std::size_t num_samples = 0;

        for (ProjectedEdgeInfo const & projected_edge_info : projected_edge_infos) {
            TexturePatch const & texture_patch = texture_patches[projected_edge_info.texture_patch_id];
            const int texture_patch_label = texture_patch.get_label();
            if (texture_patch_label == label1 || texture_patch_label == label2){
                if (texture_patch_label == label1)
                    color1_accum.add(sample_edge(texture_patch, projected_edge_info.p1, projected_edge_info.p2), length);

                if (texture_patch_label == label2)
                    color2_accum.add(sample_edge(texture_patch, projected_edge_info.p1, projected_edge_info.p2), length);

                num_samples++;
            }
        }
        assert(num_samples == 2);
    }

    math::Vec3f color1 = color1_accum.normalized();
    math::Vec3f color2 = color2_accum.normalized();

    /* The order is essential. */
    math::Vec3f difference = color2 - color1;

    assert(!isnan(difference[0]) && !isnan(difference[1]) && !isnan(difference[2]));

    return difference;
}

void
global_seam_leveling(UniGraph const & graph, mve::TriangleMesh::ConstPtr mesh,
    mve::VertexInfoList::ConstPtr vertex_infos,
    std::vector<std::vector<VertexProjectionInfo> > const & vertex_projection_infos,
    std::vector<TexturePatch> * texture_patches) {

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

        std::vector<std::size_t> faces = vertex_infos->at(i).faces;
        std::set<std::size_t>::iterator it = label_set.begin();
        for (std::size_t j = 0; j < faces.size(); ++j) {
            std::size_t label = graph.get_label(faces[j]);
            label_set.insert(it, label);
        }

        for (it = label_set.begin(); it != label_set.end(); ++it) {
            std::size_t label = *it;
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
            std::vector<std::size_t> const & adj_verts = vertex_infos->at(i).verts;
            for (std::size_t k = 0; k < adj_verts.size(); ++k) {
                std::size_t adj_vertex = adj_verts[k];
                for (std::size_t l = 0; l < labels[adj_vertex].size(); ++l) {
                    std::size_t label = labels[i][j];
                    std::size_t adj_vertex_label = labels[adj_vertex][l];
                    if (label != 0 && i < adj_vertex && label == adj_vertex_label) {
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
                if (label1 != 0 && label2 != 0 && label1 < label2) {

                    std::vector<MeshEdge> seam_edges;
                    find_seam_edges_for_vertex_label_combination(graph, mesh, vertex_infos, i, label1, label2, &seam_edges);

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

    SpMat Lhs = A.transpose() * A + Gamma.transpose() * Gamma;
    /* Only keep lower triangle (CG only uses the lower), prune the rest and compress matrix. */
    Lhs.prune([](const int& row, const int& col, const float& value) -> bool {
            return col <= row && value != 0.0f;
        }); // value != 0.0f is only to suppress a compiler warning

    std::vector<std::map<std::size_t, math::Vec3f> > adjust_values(num_vertices);
    std::cout << " done." << std::endl;
    std::cout << "\tLhs dimensionality: " << Lhs.rows() << " x " << Lhs.cols() << std::endl;

    util::WallTimer timer;
    std::cout << "Calculating adjustments:"<< std::endl;
    #pragma omp parallel for
    for (std::size_t channel = 0; channel < 3; ++channel) {
        /* Prepare solver. */
        Eigen::ConjugateGradient<SpMat, Eigen::Lower> cg;
        cg.setMaxIterations(1000);
        cg.setTolerance(0.0001);
        cg.compute(Lhs);

        /* Prepare right hand side. */
        Eigen::VectorXf b(A_rows);
        for (std::size_t i = 0; i < coefficients_b.size(); ++i)
            b[i] = coefficients_b[i][channel];
        Eigen::VectorXf Rhs = SpMat(A.transpose()) * b;

        /* Solve for x. */
        Eigen::VectorXf x(x_rows);
        x = cg.solve(Rhs);
        #pragma omp critical
        std::cout << "\t\tColor channel " << channel << ": CG took "
            << cg.iterations() << " iterations. Residual is " << cg.error() << std::endl;

        /* Subtract mean because system is underconstrained and we seek the solution with minimal adjustments. */
        x = x.array() - x.mean();

        #pragma omp critical
        {
            for (std::size_t i = 0; i < num_vertices; ++i) {
                for (std::size_t j = 0; j < labels[i].size(); ++j) {
                    std::size_t label = labels[i][j];
                    adjust_values[i][label][channel] = x[vertlabel2row[i][label]];
                }
            }
        }
    }
    std::cout << "\t\tTook " << timer.get_elapsed_sec() << " seconds" << std::endl;


    mve::TriangleMesh::FaceList const & mesh_faces = mesh->get_faces();

    ProgressCounter texture_patch_counter("\tAdjusting texture patches", texture_patches->size());
    #pragma omp parallel for schedule(dynamic)
    for (std::size_t i = 0; i < texture_patches->size(); ++i) {
        texture_patch_counter.progress<SIMPLE>();

        TexturePatch * texture_patch = &texture_patches->at(i);

        int label = texture_patch->get_label();
        std::vector<std::size_t> const & faces = texture_patch->get_faces();

        std::vector<math::Vec3f> patch_adjust_values(faces.size() * 3);
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

math::Vec3f
mean_color_of_edge_point(std::vector<ProjectedEdgeInfo> projected_edge_infos,
    std::vector<TexturePatch> const & texture_patches, float t) {

    assert(0.0f <= t && t <= 1.0f);
    math::Accum<math::Vec3f> color_accum(math::Vec3f(0.0f));

    for (ProjectedEdgeInfo const & projected_edge_info : projected_edge_infos) {
        math::Vec2f pixel = projected_edge_info.p1 * t + (1.0f - t) * projected_edge_info.p2;
        math::Vec3f color = texture_patches[projected_edge_info.texture_patch_id].get_pixel_value(pixel);
        color_accum.add(color, 1.0f);
    }

    math::Vec3f mean_color = color_accum.normalized();
    return mean_color;
}

void
draw_line(math::Vec2f p1, math::Vec2f p2, std::vector<ProjectedEdgeInfo> const & projected_edge_infos,
    std::vector<TexturePatch> const & texture_patches, TexturePatch * texture_patch) {
    /* http://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm */

    int x0 = std::floor(p1[0] + 0.5f);
    int y0 = std::floor(p1[1] + 0.5f);
    int const x1 = std::floor(p2[0] + 0.5f);
    int const y1 = std::floor(p2[1] + 0.5f);

    float tdx = static_cast<float>(x1 - x0);
    float tdy = static_cast<float>(y1 - y0);
    float length = std::sqrt(tdx * tdx + tdy * tdy);

    int const dx = std::abs(x1 - x0);
    int const dy = std::abs(y1 - y0) ;
    int const sx = x0 < x1 ? 1 : -1;
    int const sy = y0 < y1 ? 1 : -1;
    int err = dx - dy;

    int x = x0;
    int y = y0;
    while (true) {
        math::Vec2i pixel(x, y);

        tdx = static_cast<float>(x1 - x);
        tdy = static_cast<float>(y1 - y);

        /* If the length is zero we sample the midpoint of the projected edge. */
        float t = length != 0.0f ? std::sqrt(tdx * tdx + tdy * tdy) / length : 0.5f;

        texture_patch->set_pixel_value(pixel, mean_color_of_edge_point(projected_edge_infos, texture_patches, t));
        if (x == x1 && y == y1)
            break;

        int const e2 = 2 * err;
        if (e2 > -dy) {
            err -= dy;
            x += sx;
        }
        if (e2 < dx) {
            err += dx;
            y += sy;
        }
    }
}


void
local_seam_leveling(UniGraph const & graph, mve::TriangleMesh::ConstPtr mesh,
    std::vector<std::vector<VertexProjectionInfo> > const & vertex_projection_infos,
    std::vector<TexturePatch> * texture_patches) {

    std::vector<TexturePatch> orig_texture_patches(texture_patches->begin(), texture_patches->end());

    std::size_t const num_vertices = vertex_projection_infos.size();

    {
        std::vector<MeshEdge> seam_edges;
        find_seam_edges(graph, mesh, &seam_edges);

        for (MeshEdge seam_edge : seam_edges){
            std::vector<ProjectedEdgeInfo> projected_edge_infos;
            find_mesh_edge_projections(vertex_projection_infos, seam_edge, &projected_edge_infos);

            for (ProjectedEdgeInfo const & projected_edge_info : projected_edge_infos) {
                draw_line(projected_edge_info.p1, projected_edge_info.p2, projected_edge_infos,
                    orig_texture_patches, &(texture_patches->at(projected_edge_info.texture_patch_id)));
            }
        }
    }

    for (std::size_t i = 0; i < num_vertices; ++i) {
        std::vector<VertexProjectionInfo> const & projection_infos = vertex_projection_infos[i];
        if (projection_infos.size() <= 1) continue;

        math::Accum<math::Vec3f> color_accum(math::Vec3f(0.0f));
        for (std::size_t j = 0; j < projection_infos.size(); ++j) {
            VertexProjectionInfo const & projection_info = projection_infos[j];
            TexturePatch const & original_texture_patch = orig_texture_patches[projection_info.texture_patch_id];
            math::Vec3f color = original_texture_patch.get_pixel_value(projection_info.projection);
            color_accum.add(color, 1.0f);
        }

        math::Vec3f mean_color = color_accum.normalized();

        for (std::size_t j = 0; j < projection_infos.size(); ++j){
            VertexProjectionInfo const & projection_info = projection_infos[j];
            math::Vec2i pixel(projection_info.projection +  math::Vec2f(0.5f, 0.5f));
            TexturePatch * texture_patch = &(texture_patches->at(projection_info.texture_patch_id));
            texture_patch->set_pixel_value(pixel, mean_color);
        }
    }

    ProgressCounter texture_patch_counter("\tBlending texture patches", texture_patches->size());
    #pragma omp parallel for schedule(dynamic)
    for (std::size_t i = 0; i < texture_patches->size(); ++i) {
        TexturePatch * texture_patch = &texture_patches->at(i);
        texture_patch_counter.progress<SIMPLE>();
        texture_patch->prepare_blending_mask(STRIP_SIZE);
        texture_patch->blend(orig_texture_patches[i].get_image());
        texture_patch->release_blending_mask();
        texture_patch->erode_validity_mask();
        texture_patch_counter.inc();
    }
}
