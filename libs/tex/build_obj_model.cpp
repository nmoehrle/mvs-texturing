/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <iostream>

#include "defines.h"
#include "texture_atlas.h"
#include "obj_model.h"

TEX_NAMESPACE_BEGIN

void
build_model(mve::TriangleMesh::ConstPtr mesh,
    std::vector<TextureAtlas::Ptr> const & texture_atlases, ObjModel * obj_model)  {

    mve::TriangleMesh::VertexList const & mesh_vertices = mesh->get_vertices();
    mve::TriangleMesh::NormalList const & mesh_normals = mesh->get_vertex_normals();
    mve::TriangleMesh::FaceList const & mesh_faces = mesh->get_faces();

    ObjModel::Vertices & vertices = obj_model->get_vertices();
    vertices.insert(vertices.begin(), mesh_vertices.begin(), mesh_vertices.end());
    ObjModel::Normals & normals = obj_model->get_normals();
    normals.insert(normals.begin(), mesh_normals.begin(), mesh_normals.end());
    ObjModel::TexCoords & texcoords = obj_model->get_texcoords();

    ObjModel::Groups & groups = obj_model->get_groups();
    MaterialLib & material_lib = obj_model->get_material_lib();

    for (TextureAtlas::Ptr texture_atlas : texture_atlases) {

        Material material;
        const std::size_t n = material_lib.size();
        material.name = std::string("material") + util::string::get_filled(n, 4);
        material.diffuse_map = texture_atlas->get_image();
        material_lib.push_back(material);

        groups.push_back(ObjModel::Group());
        ObjModel::Group & group = groups.back();
        group.material_name = material.name;

        TextureAtlas::Faces const & atlas_faces = texture_atlas->get_faces();
        TextureAtlas::Texcoords const & atlas_texcoords = texture_atlas->get_texcoords();
        TextureAtlas::TexcoordIds const & atlas_texcoord_ids = texture_atlas->get_texcoord_ids();

        std::size_t texcoord_id_offset = texcoords.size();

        texcoords.insert(texcoords.end(), atlas_texcoords.begin(),
            atlas_texcoords.end());

        for (std::size_t i = 0; i < atlas_faces.size(); ++i) {
            std::size_t mesh_face_pos = atlas_faces[i] * 3;

            std::size_t vertex_ids[] = {
                mesh_faces[mesh_face_pos],
                mesh_faces[mesh_face_pos + 1],
                mesh_faces[mesh_face_pos + 2]
            };
            std::size_t * normal_ids = vertex_ids;

            std::size_t texcoord_ids[] = {
                texcoord_id_offset + atlas_texcoord_ids[i * 3],
                texcoord_id_offset + atlas_texcoord_ids[i * 3 + 1],
                texcoord_id_offset + atlas_texcoord_ids[i * 3 + 2]
            };

            group.faces.push_back(ObjModel::Face());
            ObjModel::Face & face = group.faces.back();
            std::copy(vertex_ids, vertex_ids + 3, face.vertex_ids);
            std::copy(texcoord_ids, texcoord_ids + 3, face.texcoord_ids);
            std::copy(normal_ids, normal_ids + 3, face.normal_ids);
        }
    }
    //TODO remove unreferenced vertices/normals.
}

TEX_NAMESPACE_END
