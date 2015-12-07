/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <set>

#include "seam_leveling.h"

TEX_NAMESPACE_BEGIN

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
            if (label1 == label2) continue;

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
find_mesh_edge_projections(
    std::vector<std::vector<VertexProjectionInfo> > const & vertex_projection_infos,
    MeshEdge mesh_edge, std::vector<EdgeProjectionInfo> * edge_projection_infos) {
    std::vector<VertexProjectionInfo> const & v1_projection_infos = vertex_projection_infos[mesh_edge.v1];
    std::vector<VertexProjectionInfo> const & v2_projection_infos = vertex_projection_infos[mesh_edge.v2];

    /* Use a set to eliminate duplicates which may occur if the mesh is degenerated. */
    std::set<EdgeProjectionInfo> edge_projection_infos_set;

    for (VertexProjectionInfo v1_projection_info : v1_projection_infos) {
        for (VertexProjectionInfo v2_projection_info : v2_projection_infos) {
            if (v1_projection_info.texture_patch_id != v2_projection_info.texture_patch_id) continue;

            for (std::size_t face_id1 : v1_projection_info.faces) {
                for (std::size_t face_id2 : v2_projection_info.faces) {
                    if(face_id1 != face_id2) continue;

                    std::size_t texture_patch_id = v1_projection_info.texture_patch_id;
                    math::Vec2f p1 = v1_projection_info.projection;
                    math::Vec2f p2 = v2_projection_info.projection;

                    EdgeProjectionInfo edge_projection_info = {texture_patch_id, p1, p2};
                    edge_projection_infos_set.insert(edge_projection_info);
                }
            }
        }
    }

    edge_projection_infos->insert(edge_projection_infos->end(), edge_projection_infos_set.begin(), edge_projection_infos_set.end());
}

TEX_NAMESPACE_END
