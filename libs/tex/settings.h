/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef TEX_SETTINGS_HEADER
#define TEX_SETTINGS_HEADER

#include <sstream>
#include <vector>
#include <string>

#include "defines.h"

template <typename T>
const std::vector<std::string> choice_strings();

template <typename T> inline
const std::string choice_string(T i) {
    return choice_strings<T>()[static_cast<std::size_t>(i)];
}

template <typename T> inline
const std::string choices() {
    const std::vector<std::string> strings = choice_strings<T>();
    const std::size_t n = strings.size();
    std::stringstream ss;
    for (std::size_t i = 0; i < n; ++i) {
        ss << strings[i];
        if(i != n - 1) ss << ", ";
    }
    return ss.str();
}

template <typename T> inline
T parse_choice(std::string s) {
    const std::vector<std::string> strings = choice_strings<T>();
    const std::size_t n = strings.size();
    for (std::size_t i = 0; i < n; ++i) {
        if (s == strings[i]) {
            return static_cast<T>(i);
        }
    }

    std::stringstream ss;
    ss << "Invalid choice: " << s << " (Available choices: ";
    ss << choices<T>() << ")";

    throw std::invalid_argument(ss.str());
}

TEX_NAMESPACE_BEGIN

/** Enum representing a data term. */
enum DataTerm {
    AREA = 0,
    GMI = 1
};

/** Enum representing a smoothness term. */
enum SmoothnessTerm {
    POTTS = 0
};

/** Enum representing the choice of outlier removal. */
enum OutlierRemoval {
    NONE = 0,
    GAUSS_DAMPING = 1,
    GAUSS_CLAMPING = 2
};

struct Settings {
    bool verbose;

    DataTerm data_term;
    SmoothnessTerm smoothness_term;
    OutlierRemoval outlier_removal;

    bool geometric_visibility_test;
    bool global_seam_leveling;
    bool local_seam_leveling;
    bool hole_filling;
    bool keep_unseen_faces;
};

TEX_NAMESPACE_END

template <> inline
const std::vector<std::string> choice_strings<tex::DataTerm>() {
    return {"area", "gmi"};
}

template <> inline
const std::vector<std::string> choice_strings<tex::SmoothnessTerm>() {
    return {"potts"};
}

template <> inline
const std::vector<std::string> choice_strings<tex::OutlierRemoval>() {
    return {"none", "gauss_damping", "gauss_clamping"};
}

#endif /* TEX_SETTINGS_HEADER */
