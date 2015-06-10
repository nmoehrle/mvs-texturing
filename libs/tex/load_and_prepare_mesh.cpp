#include "texturing.h"
#include "mve/mesh_io_ply.h"

mve::TriangleMesh::Ptr
load_and_prepare_mesh(const std::string & filename) {
    mve::TriangleMesh::Ptr mesh;
    try {
        mesh = mve::geom::load_ply_mesh(filename);
    } catch (std::exception& e) {
        std::cerr << "\tCould not load mesh: "<< e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    mve::VertexInfoList::Ptr vertex_infos = mve::VertexInfoList::create(mesh);

    /* Remove redundant faces. */
    int num_redundant = 0;
    mve::TriangleMesh::FaceList & faces = mesh->get_faces();
    mve::TriangleMesh::FaceList new_faces;
    new_faces.reserve(faces.size());
    for (std::size_t i = 0; i < faces.size(); i += 3) {
        std::size_t face_id = i / 3;
        bool redundant = false;
        for (std::size_t j = 0; !redundant && j < 3; ++j) {
            mve::MeshVertexInfo::FaceRefList const & adj_faces = vertex_infos->at(faces[i + j]).faces;
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

        if (redundant)
            ++num_redundant;
        else
            new_faces.insert(new_faces.end(), &faces[i], &faces[i + 3]);
    }

    faces.swap(new_faces);

    if (num_redundant > 0)
        std::cout << "\tRemoved " << num_redundant << " redundant faces." << std::endl;

    /* Ensure face and vertex normals. */
    mesh->ensure_normals(true, true);

    return mesh;
}
