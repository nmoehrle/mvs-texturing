/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include "texturing.h"

#include "progress_counter.h"

TEX_NAMESPACE_BEGIN

void
build_adjacency_graph(mve::TriangleMesh::ConstPtr mesh,
    mve::MeshInfo const & mesh_info, UniGraph * graph)  {

    mve::TriangleMesh::FaceList const & faces = mesh->get_faces();
    std::size_t const num_faces = faces.size() / 3;

    ProgressCounter face_counter("\tAdding edges", num_faces);
    for (std::size_t i = 0; i < faces.size(); i += 3) {
        face_counter.progress<SIMPLE>();

        std::size_t v1 = faces[i];
        std::size_t v2 = faces[i + 1];
        std::size_t v3 = faces[i + 2];

        std::vector<std::size_t> adj_faces;
        mesh_info.get_faces_for_edge(v1, v2, &adj_faces);
        mesh_info.get_faces_for_edge(v2, v3, &adj_faces);
        mesh_info.get_faces_for_edge(v3, v1, &adj_faces);

        for (std::size_t j = 0; j < adj_faces.size(); ++j) {
            /* Face id vs. face position. */
            std::size_t face = i / 3;
            std::size_t adj_face = adj_faces[j];

            /* Avoid self referencing. */
            if (face != adj_face) {
                /* Edge not already in graph? */
                if (!graph->has_edge(face, adj_face)){
                    graph->add_edge(face, adj_face);
                }
            }
        }
        face_counter.inc();
    }

    std::cout << "\t" << graph->num_edges() << " total edges." << std::endl;
}

TEX_NAMESPACE_END
