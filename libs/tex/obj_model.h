/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef TEX_OBJMODEL_HEADER
#define TEX_OBJMODEL_HEADER

#include "material_lib.h"

/**
  * Class representing a obj model.
  */
class ObjModel {
public:
    struct Face {
        std::size_t vertex_ids[3];
        std::size_t texcoord_ids[3];
        std::size_t normal_ids[3];
    };

    struct Group {
        std::string material_name;
        std::vector<Face> faces;
    };

    typedef std::vector<math::Vec3f> Vertices;
    typedef std::vector<math::Vec2f> TexCoords;
    typedef std::vector<math::Vec3f> Normals;
    typedef std::vector<Group> Groups;

private:
    Vertices vertices;
    TexCoords texcoords;
    Normals normals;
    Groups groups;
    MaterialLib material_lib;

public:
    /** Saves the obj model to an .obj file, its material lib and the materials with the given prefix. */
    void save_to_files(std::string const & prefix) const;
    ObjModel();

    MaterialLib & get_material_lib(void);
    Vertices & get_vertices(void);
    TexCoords & get_texcoords(void);
    Normals & get_normals(void);
    Groups & get_groups(void);

    static void save(ObjModel const & model, std::string const & prefix);
};

inline
MaterialLib &
ObjModel::get_material_lib(void) {
    return material_lib;
}

inline
ObjModel::Vertices &
ObjModel::get_vertices(void) {
    return vertices;
}

inline
ObjModel::TexCoords &
ObjModel::get_texcoords(void) {
    return texcoords;
}

inline
ObjModel::Normals &
ObjModel::get_normals(void) {
    return normals;
}

inline
ObjModel::Groups &
ObjModel::get_groups(void) {
    return groups;
}

#endif /* TEX_OBJMODEL_HEADER */
