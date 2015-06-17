#pragma once

#include "mve/image.h"

/**
  * Class representing a lambertian material.
  */
class Material {
    private:
        mve::ByteImage::Ptr diffuse_map;
        //TODO: implement BRDF properties

    public:
        Material(mve::ByteImage::Ptr diffuse_map);

        mve::ByteImage::Ptr get_diffuse_map() const;
};

inline
mve::ByteImage::Ptr
Material::get_diffuse_map() const {
    return diffuse_map;
}
