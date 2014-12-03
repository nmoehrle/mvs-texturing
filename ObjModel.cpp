#include "ObjModel.h"

ObjModel::ObjModel() {
}

void
ObjModel::save_to_files(std::string const & prefix) {
    material_lib.save_to_files(prefix);

    std::string name = util::fs::basename(prefix);
    std::ofstream out((prefix + ".obj").c_str());
    if (!out.good())
        throw util::FileException(prefix + ".obj", std::strerror(errno));

    out << "mtllib " << name << ".mtl" << std::endl;

    out << std::fixed << std::setprecision(6);
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        out << "v " << vertices[i][0] << " "
            << vertices[i][1] << " "
            << vertices[i][2] << std::endl;
    }

    for (std::size_t i = 0; i < texcoords.size(); ++i) {
        out << "vt " << texcoords[i][0] << " "
            << texcoords[i][1] << std::endl;
    }

    for (std::size_t i = 0; i < normals.size(); ++i) {
        out << "vn " << normals[i][0] << " "
            << normals[i][1] << " "
            << normals[i][2] << std::endl;
    }

    for (std::size_t i = 0; i < groups.size(); ++i) {
        out << "usemtl " << groups[i].material_name << std::endl;
        for (std::size_t j = 0; j < groups[i].faces.size(); ++j) {
            Face & face =  groups[i].faces[j];
            out << "f";
            for (std::size_t k = 0; k < 3; ++k) {
                out << " " << face.vertices[k]  + OBJ_INDEX_OFFSET
                    << "/" << face.texcoords[k]  + OBJ_INDEX_OFFSET
                    << "/" << face.normals[k]  + OBJ_INDEX_OFFSET;
            }
            out << std::endl;
        }
    }
    out.close();
}
