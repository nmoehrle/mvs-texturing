/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef TEX_TEXTUREATLAS_HEADER
#define TEX_TEXTUREATLAS_HEADER


#include <vector>

#include <util/exception.h>
#include <math/vector.h>
#include <mve/mesh.h>
#include <mve/image.h>

#include "tri.h"
#include "texture_patch.h"
#include "rectangular_bin.h"

/**
  * Class representing a texture atlas.
  */
class TextureAtlas {
    public:
        typedef std::shared_ptr<TextureAtlas> Ptr;

        typedef std::vector<std::size_t> Faces;
        typedef std::vector<std::size_t> TexcoordIds;
        typedef std::vector<math::Vec2f> Texcoords;

    private:
        unsigned int const size;
        unsigned int const padding;
        bool finalized;

        Faces faces;
        Texcoords texcoords;
        TexcoordIds texcoord_ids;

        mve::ImageBase::Ptr image;
        mve::ByteImage::Ptr validity_mask;

        RectangularBin::Ptr bin;

        template <typename T>
        void apply_edge_padding(void);

        void merge_texcoords(void);

    public:
        TextureAtlas(unsigned int size, mve::ImageType type);

        static TextureAtlas::Ptr create(unsigned int size, mve::ImageType type);

        Faces const & get_faces(void) const;
        TexcoordIds const & get_texcoord_ids(void) const;
        Texcoords const & get_texcoords(void) const;
        mve::ImageBase::Ptr get_image(void) const;

        bool insert(TexturePatch::ConstPtr texture_patch);

        void finalize(void);
};

inline TextureAtlas::Ptr
TextureAtlas::create(unsigned int size, mve::ImageType type) {
    return Ptr(new TextureAtlas(size, type));
}

inline TextureAtlas::Faces const &
TextureAtlas::get_faces(void) const {
    return faces;
}

inline TextureAtlas::TexcoordIds const &
TextureAtlas::get_texcoord_ids(void) const {
    return texcoord_ids;
}

inline TextureAtlas::Texcoords const &
TextureAtlas::get_texcoords(void) const {
    return texcoords;
}

inline mve::ImageBase::Ptr
TextureAtlas::get_image(void) const {
    if (!finalized) {
        throw util::Exception("Texture atlas not finalized");
    }
    return image;
}

#endif /* TEX_TEXTUREATLAS_HEADER */
