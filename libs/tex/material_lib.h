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

struct Material {
    std::string diffuse_map;
};

/**
  * Class representing a material lib of and obj model.
  */
class MaterialLib {
    private:
        std::vector<Material> materials;
        std::vector<std::string> material_names;
    public:
        MaterialLib();

        void add_material(std::string const & name, Material material);
        std::size_t size();

        /** Saves the material lib to an .mtl file and all textures of its
          * materials with the given prefix.
          */
        void save_to_files(std::string const & prefix) const;
};

inline std::size_t
MaterialLib::size() {
    return materials.size();
}

#endif /* TEX_MATERIALLIB_HEADER */
