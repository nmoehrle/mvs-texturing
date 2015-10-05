/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef TEX_MATERIALLIB_HEADER
#define TEX_MATERIALLIB_HEADER

#include <vector>
#include <mve/image.h>

struct Material {
    std::string name;
    mve::ByteImage::ConstPtr diffuse_map;
};

/**
  * Class representing a material lib of and obj model.
  */
class MaterialLib : public std::vector<Material>{

    public:
        /** Saves the material lib to an .mtl file and all maps of its
          * materials with the given prefix.
          */
        void save_to_files(std::string const & prefix) const;
};

#endif /* TEX_MATERIALLIB_HEADER */
