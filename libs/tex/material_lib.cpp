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
#include <mve/image_io.h>

#include "material_lib.h"

void
MaterialLib::save_to_files(std::string const & prefix) const {
    std::string filename = prefix + ".mtl";
    std::ofstream out(filename.c_str());
    if (!out.good())
        throw util::FileException(filename, std::strerror(errno));

    std::string const name = util::fs::basename(prefix);

    for (Material const & material : *this) {
        std::string diffuse_map_postfix = "_" + material.name + "_map_Kd.png";
        out << "newmtl " << material.name << std::endl
            << "Ka 1.000000 1.000000 1.000000" << std::endl
            << "Kd 1.000000 1.000000 1.000000" << std::endl
            << "Ks 0.000000 0.000000 0.000000" << std::endl
            << "Tr 1.000000" << std::endl
            << "illum 1" << std::endl
            << "Ns 1.000000" << std::endl
            << "map_Kd " << name + diffuse_map_postfix << std::endl;
    }
    out.close();

    for (Material const & material : *this) {
        std::string filename = prefix + "_" + material.name + "_map_Kd.png";
        mve::image::save_png_file(material.diffuse_map, filename);
    }
}
