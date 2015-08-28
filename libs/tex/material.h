/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef TEX_MATERIAL_HEADER
#define TEX_MATERIAL_HEADER

#include <mve/image.h>

/**
  * Class representing a lambertian material.
  */
class Material {
    private:
        mve::ByteImage::ConstPtr diffuse_map;
        //TODO: implement BRDF properties

    public:
        Material(mve::ByteImage::ConstPtr diffuse_map);

        mve::ByteImage::ConstPtr get_diffuse_map() const;
};

inline
mve::ByteImage::ConstPtr
Material::get_diffuse_map() const {
    return diffuse_map;
}

inline
Material::Material(mve::ByteImage::ConstPtr diffuse_map)
    : diffuse_map(diffuse_map) {
}

#endif /* TEX_MATERIAL_HEADER */
