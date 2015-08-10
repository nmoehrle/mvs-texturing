#define MAX_TEXTURE_SIZE (8*1024)
#include "texturing.h"
#include <set>

TEX_NAMESPACE_BEGIN

/**
  * Heuristic to calculate an appropriate texture atlas size.
  * @warning asserts that no texture patch exceeds the dimensions
  * of the maximal possible texture atlas size.
  */
int
calculate_texture_size(std::list<TexturePatch> & texture_patches) {
    double area = 0;
    int max_width = 0;
    int max_height = 0;
    std::list<TexturePatch>::iterator it = texture_patches.begin();
    for (; it != texture_patches.end(); it++) {
        int const width = it->get_width();
        int const height = it->get_height();
        area += width * height;
        max_width = std::max(max_width, width);
        max_height = std::max(max_height, height);
    }

    /* A TexturePatch must fit into the largest possible texture. */
    int const max_padding = MAX_TEXTURE_SIZE >> 7;
    int const max_patch_size = MAX_TEXTURE_SIZE - 2 * max_padding;
    (void) max_patch_size; //Suppress unused variable if -DNDEBUG
    assert(max_patch_size >= max_width && max_patch_size >= max_height);

    /* Heuristic to determine a proper texture size. */
    /* Approximate area which will be needed */
    double approx_area = 1.4 * area;

    /* Probably fits into a single texture atlas. */
    if (approx_area < MAX_TEXTURE_SIZE * MAX_TEXTURE_SIZE) {
        /* Next larger power of two. */
        int size = pow(2, ceil(log(sqrt(approx_area)) / log(2.0)));

        while (size < max_width || size < max_height) {
            size *= 2;
        }

        return size;
    } else {
        return MAX_TEXTURE_SIZE;
    }
}

/**
  * Copies the src image into the dest image at the given position,
  * optionally adding a border.
  * @warning asserts that the given src image fits into the given dest image.
  */
void copy_into(mve::ByteImage::ConstPtr src, int x, int y,
    mve::ByteImage::Ptr dest, int border = 0) {

    assert(x >= 0 && x + src->width() + 2 * border <= dest->width());
    assert(y >= 0 && y + src->height() + 2 * border <= dest->height());

    for (int i = 0; i < src->width() + 2 * border; ++i) {
        for(int j = 0; j < src->height() + 2 * border; j++) {
            int sx = i - border;
            int sy = j - border;

            if (sx < 0 || sx >= src->width() || sy < 0 || sy >= src->height())
                continue;

            for (int c = 0; c < src->channels(); ++c)
                dest->at(x + i, y + j, c) = src->at(sx, sy, c);
        }
    }
}

/**
  * Iteratively dilates all valid pixels using a 3x3 gaussian kernel,
  * to determine the value for the new pixels.
  * @warning asserts matching dimensions of image and validity_mask
  * @warning asserts a three channel image
  * @warning asserts a one channel validity_mask
  */
void
dilate_valid_pixel(mve::ByteImage::Ptr image, mve::ByteImage::Ptr validity_mask) {

    assert(image->width() == validity_mask->width());
    assert(image->height() == validity_mask->height());
    assert(image->channels() == 3);
    assert(validity_mask->channels() == 1);

    const int width = image->width();
    const int height = image->height();
    const int size = width;

    math::Matrix<float, 3, 3> gauss;
    gauss[0] = 1.0f; gauss[1] = 2.0f; gauss[2] = 1.0f;
    gauss[3] = 2.0f; gauss[4] = 4.0f; gauss[5] = 2.0f;
    gauss[6] = 1.0f; gauss[7] = 2.0f; gauss[8] = 1.0f;
    gauss /= 16.0f;

    /* Calculate the set of invalid pixels at the border of texture patches. */
    PixelSet invalid_border_pixels;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            if (validity_mask->at(x, y, 0) == 255) continue;

            /* Check the direct neighbourhood of all invalid pixels. */
            bool at_border = false;
            for (int j = -1; j <= 1 && !at_border; ++j) {
                for (int i = -1; i <= 1 && !at_border; ++i) {
                    int nx = x + i;
                    int ny = y + j;
                    /* If the invalid pixel has a valid neighbour: */
                    if (0 <= nx && nx < width &&
                        0 <= ny && ny < height &&
                        validity_mask->at(nx, ny, 0) == 255) {

                        /* Add the pixel to the set of invalid border pixels. */
                        invalid_border_pixels.insert(std::pair<int, int>(x, y));
                        at_border = true;
                    }
                }
            }
        }
    }

    mve::ByteImage::Ptr new_validity_mask = validity_mask->duplicate();

    /* Iteratively dilate border pixels until padding constants are reached. */
    for (int n = 0; n < (size >> 7); ++n) {
        PixelVector new_valid_pixels;

        PixelSet::iterator it = invalid_border_pixels.begin();
        for (;it != invalid_border_pixels.end(); it++) {
            int x = it->first;
            int y = it->second;

            bool now_valid = false;
            /* Calculate new pixel value. */
            for (int c = 0; c < 3; ++c) {
                float norm = 0.0f;
                float value = 0.0f;
                for (int j = -1; j <= 1; ++j) {
                    for (int i = -1; i <= 1; ++i) {
                        int nx = x + i;
                        int ny = y + j;
                        if (0 <= nx && nx < width &&
                            0 <= ny && ny < height &&
                            new_validity_mask->at(nx, ny, 0) == 255) {

                            float w = gauss[(j + 1) * 3 + (i + 1)];
                            norm += w;
                            value += (image->at(nx, ny, c) / 255.0f) * w;
                        }
                    }
                }

                if (norm == 0.0f)
                    continue;

                now_valid = true;
                image->at(x, y, c) = (value / norm) * 255.0f;
            }

            if (now_valid) {
                new_valid_pixels.push_back(*it);
            }
        }

        invalid_border_pixels.clear();

        /* Mark the new valid pixels valid in the validity mask. */
        for (std::size_t i = 0; i < new_valid_pixels.size(); ++i) {
             int x = new_valid_pixels[i].first;
             int y = new_valid_pixels[i].second;

             new_validity_mask->at(x, y, 0) = 255;
        }

        /* Calculate the set of invalid pixels at the border of the valid area. */
        for (std::size_t i = 0; i < new_valid_pixels.size(); ++i) {
            int x = new_valid_pixels[i].first;
            int y = new_valid_pixels[i].second;

            for (int j = -1; j <= 1; ++j) {
                 for (int i = -1; i <= 1; ++i) {
                     int nx = x + i;
                     int ny = y + j;
                     if (0 <= nx && nx < width &&
                         0 <= ny && ny < height &&
                         new_validity_mask->at(nx, ny, 0) == 0) {

                         invalid_border_pixels.insert(std::pair<int, int>(nx, ny));
                    }
                }
            }
        }
    }
}


void
build_model(mve::TriangleMesh::ConstPtr mesh,
    std::vector<TexturePatch> const & _texture_patches, ObjModel * obj_model)  {

    mve::TriangleMesh::VertexList const & mesh_vertices = mesh->get_vertices();
    mve::TriangleMesh::NormalList const & mesh_normals = mesh->get_vertex_normals();
    mve::TriangleMesh::FaceList const & mesh_faces = mesh->get_faces();
    std::size_t num_faces = mesh_faces.size() / 3;

    ObjModel::Vertices & vertices = obj_model->get_vertices();
    vertices.insert(vertices.begin(), mesh_vertices.begin(), mesh_vertices.end());
    ObjModel::Normals & normals = obj_model->get_normals();
    normals.insert(normals.begin(), mesh_normals.begin(), mesh_normals.end());
    ObjModel::TexCoords & texcoords = obj_model->get_texcoords();
    /* Preallocate for the maximum number of faces. */
    texcoords.resize(num_faces * 3);

    ObjModel::Groups & groups = obj_model->get_groups();
    MaterialLib & material_lib = obj_model->get_material_lib();

    std::size_t num_new_faces = 0;

    std::list<TexturePatch> texture_patches(_texture_patches.begin(), _texture_patches.end());

    std::cout << "\tSorting texture patches... " << std::flush;
    /* Improve the bin-packing algorithm efficiency by sorting texture patches
     * in descending order of size. */
    texture_patches.sort();
    texture_patches.reverse();
    std::cout << "done." << std::endl;

    std::cout << "\tGenerating texture atlases... " << std::flush;
    std::size_t const total_num_patches = texture_patches.size();
    std::size_t remaining_patches = texture_patches.size();
    std::ofstream tty("/dev/tty", std::ios_base::out);

    #pragma omp parallel
    {

    #pragma omp single nowait
    {

    while (!texture_patches.empty()) {
        ObjModel::Group group;

        const std::size_t n = material_lib.size();
        group.material_name = std::string("material") + util::string::get_filled(n, 4);

        int texture_size = calculate_texture_size(texture_patches);
        int padding = texture_size >> 7;
        RectangularBin bin(texture_size, texture_size);
        mve::ByteImage::Ptr texture = mve::ByteImage::create(texture_size, texture_size, 3);
        mve::ByteImage::Ptr validity_mask = mve::ByteImage::create(texture_size, texture_size, 1);

        /* Try to insert each of the texture patches into the texture atlas. */
        std::list<TexturePatch>::iterator it = texture_patches.begin();
        for (; it != texture_patches.end();) {
            TexturePatch texture_patch = *it;

            /* Progress output */
            std::size_t done_patches = total_num_patches - remaining_patches;
            if (total_num_patches > 100 && done_patches % (total_num_patches / 100) == 0) {
                tty << "\r\tGenerating texture atlases (working on atlas " << n + 1 << ") "
                    << floor((static_cast<float>(done_patches) / total_num_patches) * 100.0f + 0.5f)
                    << "%... " << std::flush;
            }

            int const width = texture_patch.get_width() + 2 * padding;
            int const height = texture_patch.get_height() + 2 * padding;
            Rect<int> rect(0, 0, width, height);
            if (bin.insert(&rect)) {
                /* Update texture atlas and its validity mask. */
                mve::ByteImage::ConstPtr patch_image = texture_patch.get_image();

                copy_into(patch_image, rect.min_x, rect.min_y, texture, padding);
                mve::ByteImage::ConstPtr patch_validity_mask = texture_patch.get_validity_mask();
                copy_into(patch_validity_mask, rect.min_x, rect.min_y, validity_mask, padding);

                TexturePatch::Faces & patch_faces = texture_patch.get_faces();
                TexturePatch::Texcoords & patch_texcoords = texture_patch.get_texcoords();

                /* Calculate the offset of the texture patches' relative texture coordinates */
                math::Vec2f offset = math::Vec2f(rect.min_x + padding, rect.min_y + padding);
                /* Calculate the final textcoords of the faces and add them to the obj model. */
                for (std::size_t i = 0; i < patch_faces.size(); ++i) {
                    std::size_t face_id = patch_faces[i];
                    std::size_t face_pos = face_id * 3;
                    std::size_t new_face_pos = num_new_faces * 3;

                    for (int j = 0; j < 3; ++j) {
                        math::Vec2f rel_texcoord(patch_texcoords[i * 3 + j]);
                        math::Vec2f texcoord = rel_texcoord + offset;

                        /* Normalize and flip Y axis */
                        texcoord[0] = texcoord[0] / (texture_size - 1);
                        texcoord[1] = 1.0f - texcoord[1] / (texture_size - 1);
                        texcoords[new_face_pos + j] = texcoord;
                    }

                    std::size_t vertex_indices[] = {
                        mesh_faces[face_pos],
                        mesh_faces[face_pos + 1],
                        mesh_faces[face_pos + 2]
                    };

                    /* Create a new ObjModel::Face */
                    std::size_t texture_indices[] = {
                        new_face_pos,
                        new_face_pos + 1,
                        new_face_pos + 2
                    };

                    std::size_t * normal_indices = vertex_indices;
                    ObjModel::Face face;
                    std::copy(vertex_indices, vertex_indices + 3, face.vertices);
                    std::copy(texture_indices, texture_indices + 3, face.texcoords);
                    std::copy(normal_indices, normal_indices + 3, face.normals);

                    group.faces.push_back(face);

                    ++num_new_faces;
                }

                it = texture_patches.erase(it);
                --remaining_patches;
            } else {
                ++it;
            }
        }

        #pragma omp task
        dilate_valid_pixel(texture, validity_mask);

        groups.push_back(group);
        material_lib.add_material(group.material_name, texture);
    }


    /* Reduce to actual size. */
    texcoords.resize(num_new_faces * 3);

    std::cout << "\r\tGenerating texture atlases (working on atlas "
        << material_lib.size() << ") 100%... done." << std::endl;
    util::WallTimer timer;
    std::cout << "\tFilling invalid pixels of the texture atlases... " << std::flush;
    #pragma omp taskwait
    std::cout << "done. (Took: " << timer.get_elapsed_sec() << "s)" << std::endl;

    /* End of parallel and single region. */
    }
    }
}

TEX_NAMESPACE_END
