/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <algorithm>
#include <cassert>
#include <fstream>
#include <cstring>
#include <cerrno>
#include <cmath>

#include <util/file_system.h>
#include <util/exception.h>

#include "histogram.h"

Histogram::Histogram(float _min, float _max, std::size_t num_bins)
    :min(_min), max(_max), num_values(0) {
    bins.resize(num_bins);
}

void
Histogram::add_value(float value) {
    float clamped_value = std::max(min, std::min(max, value));
    std::size_t index = floor(((clamped_value - min) / (max - min)) * (bins.size() - 1));
    assert(index < bins.size());
    bins[index]++;
    ++num_values;
}

void
Histogram::save_to_file(std::string const & filename) const {
    std::ofstream out(filename.c_str());
    if (!out.good())
        throw util::FileException(filename, std::strerror(errno));

    out << "Bin, Values" << std::endl;
    for (std::size_t i = 0; i < bins.size(); ++i) {
        out << i << ", " << bins[i] << std::endl;
    }
    out.close();
}

float
Histogram::get_approx_percentile(float percentile) const {
    assert(percentile >= 0.0f && percentile <= 1.0f);

    int num = 0;
    float upper_bound = min;
    for (std::size_t i = 0; i < bins.size(); ++i) {
        if (static_cast<float>(num) / num_values > percentile)
            return upper_bound;

        num += bins[i];
        upper_bound = (static_cast<float>(i) / (bins.size() - 1)) * (max - min) + min;
    }
    return max;
}
