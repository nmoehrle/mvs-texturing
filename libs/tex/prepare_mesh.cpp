/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include "texturing.h"

TEX_NAMESPACE_BEGIN

std::size_t remove_redundant_faces(mve::MeshInfo const & mesh_info, mve::TriangleMesh::Ptr mesh) {
    mve::TriangleMesh::FaceList & faces = mesh->get_faces();
    mve::TriangleMesh::FaceList new_faces;
    new_faces.reserve(faces.size());

    std::size_t num_redundant = 0;
    for (std::size_t i = 0; i < faces.size(); i += 3) {
        std::size_t face_id = i / 3;
        bool redundant = false;
        for (std::size_t j = 0; !redundant && j < 3; ++j) {
            mve::MeshInfo::AdjacentFaces const & adj_faces = mesh_info[faces[i + j]].faces;
            for (std::size_t k = 0; !redundant && k < adj_faces.size(); ++k) {
                std::size_t adj_face_id = adj_faces[k];

                /* Remove only the redundant face with smaller id. */
                if (face_id < adj_face_id) {
                    bool identical = true;
                    /* Faces are considered identical if they consist of the same vertices. */
                    for(std::size_t l = 0; l < 3; ++l) {
                        std::size_t vertex = faces[adj_face_id * 3 + l];
                        if (std::find(&faces[i], &faces[i + 3], vertex) == &faces[i + 3]) {
                            identical = false;
                            break;
                        }
                    }

                    redundant = identical;
                }
            }
        }

        if (redundant) {
            ++num_redundant;
        } else {
            new_faces.insert(new_faces.end(), faces.cbegin() + i, faces.cbegin() + i + 3);
        }
    }

    faces.swap(new_faces);

    return num_redundant;
}

void
prepare_mesh(mve::MeshInfo * mesh_info, mve::TriangleMesh::Ptr mesh) {
    std::size_t num_redundant = remove_redundant_faces(*mesh_info, mesh);
    if (num_redundant > 0) {
        std::cout << "\tRemoved " << num_redundant << " redundant faces." << std::endl;
    }

    /* Ensure face and vertex normals. */
    mesh->ensure_normals(true, true);

    /* Update vertex infos. */
    mesh_info->clear();
    mesh_info->initialize(mesh);
}

TEX_NAMESPACE_END
