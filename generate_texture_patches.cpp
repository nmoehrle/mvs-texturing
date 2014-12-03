#include "texturing.h"

void merge_vertex_projection_infos(std::vector<std::vector<VertexProjectionInfo> > * vertex_projection_infos) {
    /* Merge vertex infos within the same texture patch. */
    #pragma omp parallel for
    for (std::size_t i = 0; i < vertex_projection_infos->size(); ++i) {
        std::vector<VertexProjectionInfo> * infos = &(vertex_projection_infos->at(i));

        std::map<std::size_t, VertexProjectionInfo> info_map;
        std::map<std::size_t, VertexProjectionInfo>::iterator it;

        for (VertexProjectionInfo const & info : *infos) {
            std::size_t texture_patch_id = info.texture_patch_id;
            if((it = info_map.find(texture_patch_id)) == info_map.end())
                info_map[texture_patch_id] = info;
            else
                it->second.faces.push_back(info.faces[0]);
        }

        infos->clear();
        infos->reserve(info_map.size());
        for (it = info_map.begin(); it != info_map.end(); ++it)
            infos->push_back(it->second);
    }
}

/** Struct representing a TexturePatch candidate - final texture patches are obtained by merging candiates. */
struct TexturePatchCandidate {
    Rect<int> bounding_box;
    TexturePatch texture_patch;
};

/** Create a TexturePatchCandidate by calculating the faces' bounding box projected into the view,
  *  relative texture coordinates and extacting the texture views relevant part
  */
TexturePatchCandidate
generate_candidate(int label, TextureView const & texture_view,
    std::vector<std::size_t> const & faces, mve::TriangleMesh::ConstPtr mesh) {

    mve::ByteImage::Ptr view_image = texture_view.get_image();
    int min_x = view_image->width(), min_y = view_image->height();
    int max_x = 0, max_y = 0;

    mve::TriangleMesh::FaceList const & mesh_faces = mesh->get_faces();
    mve::TriangleMesh::VertexList const & vertices = mesh->get_vertices();

    std::vector<math::Vec2f> texcoords;
    for (std::size_t i = 0; i < faces.size(); ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            math::Vec3f vertex = vertices[mesh_faces[faces[i] * 3 + j]];
            math::Vec2f pixel = texture_view.get_pixel_coords(vertex);

            texcoords.push_back(pixel);

            min_x = std::min(static_cast<int>(std::floor(pixel[0])) - texture_patch_border, min_x);
            min_y = std::min(static_cast<int>(std::floor(pixel[1])) - texture_patch_border, min_y);
            max_x = std::max(static_cast<int>(std::ceil(pixel[0])) + texture_patch_border, max_x);
            max_y = std::max(static_cast<int>(std::ceil(pixel[1])) + texture_patch_border, max_y);
        }
    }

    /* Calculate the relative texcoords. */
    math::Vec2f min(min_x, min_y);
    for (std::size_t i = 0; i < texcoords.size(); ++i) {
        texcoords[i] = texcoords[i] - min;
    }

    int const width = max_x - min_x;
    int const height = max_y - min_y;

    mve::ByteImage::Ptr image;
    image = mve::image::crop(view_image, width, height, min_x, min_y, *math::Vec3uc(255, 0, 255));

    TexturePatchCandidate texture_patch_candidate =
        {Rect<int>(min_x, min_y, max_x, max_y), TexturePatch(label, faces, texcoords, image)};
    return texture_patch_candidate;
}

void
generate_texture_patches(UniGraph const & graph, std::vector<TextureView> const & texture_views,
    mve::TriangleMesh::ConstPtr mesh, std::vector<std::vector<VertexProjectionInfo> > * vertex_projection_infos,
    std::vector<TexturePatch> * texture_patches) {

    util::WallTimer timer;

    std::size_t const num_vertices = mesh->get_vertices().size();
    vertex_projection_infos->resize(num_vertices);

    std::size_t num_patches = 0;

    std::cout << "\tRunning... " << std::flush;
    #pragma omp parallel for schedule(dynamic)
    for (std::size_t i = 0; i < texture_views.size(); ++i) {

        std::vector<std::vector<std::size_t> > subgraphs;
        int const label = i + 1;
        graph.get_subgraphs(label, &subgraphs);

        std::list<TexturePatchCandidate> candidates;
        for (std::size_t j = 0; j < subgraphs.size(); ++j) {
            candidates.push_back(generate_candidate(label, texture_views[i], subgraphs[j], mesh));
        }

        /* Merge candidates which contain the same image content. */
        std::list<TexturePatchCandidate>::iterator it, sit;
        for (it = candidates.begin(); it != candidates.end(); ++it) {
            for (sit = candidates.begin(); sit != candidates.end();) {
                Rect<int> bounding_box = sit->bounding_box;
                if (it != sit && bounding_box.is_inside(&it->bounding_box)) {
                    TexturePatch::Faces & faces = it->texture_patch.get_faces();
                    TexturePatch::Faces & ofaces = sit->texture_patch.get_faces();
                    faces.insert(faces.end(), ofaces.begin(), ofaces.end());

                    TexturePatch::Texcoords & texcoords = it->texture_patch.get_texcoords();
                    TexturePatch::Texcoords & otexcoords = sit->texture_patch.get_texcoords();
                    math::Vec2f offset;
                    offset[0] = sit->bounding_box.min_x - it->bounding_box.min_x;
                    offset[1] = sit->bounding_box.min_y - it->bounding_box.min_y;
                    for (std::size_t i = 0; i < otexcoords.size(); ++i) {
                        texcoords.push_back(otexcoords[i] + offset);
                    }

                    sit = candidates.erase(sit);
                } else {
                    ++sit;
                }
            }
        }

        it = candidates.begin();
        for (; it != candidates.end(); ++it) {
            mve::TriangleMesh::FaceList const & mesh_faces = mesh->get_faces();
            std::size_t texture_patch_id;

            #pragma omp critical
            {
                texture_patches->push_back(it->texture_patch);
                texture_patch_id = num_patches;
                num_patches++;
            }

            std::vector<std::size_t> const & faces = it->texture_patch.get_faces();
            std::vector<math::Vec2f> const & texcoords = it->texture_patch.get_texcoords();
            for (std::size_t i = 0; i < faces.size(); ++i) {
                std::size_t const face_id = faces[i];
                std::size_t const face_pos = face_id * 3;
                for (std::size_t j = 0; j < 3; ++j) {
                    std::size_t const vertex_id = mesh_faces[face_pos  + j];
                    math::Vec2f const projection = texcoords[i * 3 + j];

                    VertexProjectionInfo info = {texture_patch_id, projection, {face_id}};

                    #pragma omp critical
                    vertex_projection_infos->at(vertex_id).push_back(info);
                }
            }
        }
    }

    merge_vertex_projection_infos(vertex_projection_infos);

    std::cout << "done. (Took " << timer.get_elapsed_sec() << "s)" << std::endl;
    std::cout << "\t" << num_patches << " texture patches." << std::endl;
}
