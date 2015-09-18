/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <fstream>
#include <cstring>
#include <cerrno>

#include <util/exception.h>
#include <util/file_system.h>
#include <mve/image.h>
#include <mve/image_io.h>

#include "material_lib.h"

MaterialLib::MaterialLib() {
}

void
MaterialLib::add_material(std::string const & name, Material material) {
    material_names.push_back(name);
    materials.push_back(material);
}

void
MaterialLib::save_to_files(std::string const & prefix) const {
    std::string filename = prefix + ".mtl";
    std::ofstream out(filename.c_str());
    if (!out.good())
        throw util::FileException(filename, std::strerror(errno));

    std::string name = util::fs::basename(prefix);

    for (std::size_t i = 0; i < materials.size(); ++i) {
        //TODO read the material parameter
        std::string diffuse_map_postfix = "_" + material_names[i] + "_map_Kd.png";
        out << "newmtl " << material_names[i] << std::endl
            << "Ka 1.000000 1.000000 1.000000" << std::endl
            << "Kd 1.000000 1.000000 1.000000" << std::endl
            << "Ks 0.000000 0.000000 0.000000" << std::endl
            << "Tr 1.000000" << std::endl
            << "illum 1" << std::endl
            << "Ns 1.000000" << std::endl
            << "map_Kd " << name + diffuse_map_postfix << std::endl;
    }
    out.close();

    for (std::size_t i = 0; i < materials.size(); ++i) {
        std::string filename = prefix + "_" + material_names[i] + "_map_Kd.png";
        util::fs::copy_file(materials[i].diffuse_map.c_str(), filename.c_str());
    }
}
