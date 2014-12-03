#pragma once

#include <sstream>

/** Enum representing the choice of outlier removal. */
enum OutlierRemoval {
    NONE = 0,
    GAUSS_DAMPING = 1,
    GAUSS_CLAMPING = 2
};

char const * const OutlierRemovalStrings [] = {"none", "gauss_damping", "gauss_clamping"};

inline
OutlierRemoval parse_outlier_removal(std::string s) {
    std::size_t n = sizeof(OutlierRemovalStrings) / sizeof(char const *);
    for (std::size_t i = 0; i < n; ++i) {
        if (s == OutlierRemovalStrings[i]) {
            return static_cast<OutlierRemoval>(i);
        }
    }

    std::stringstream ss;
    ss << "Invalid outlier removal method: " << s << " (Available outlier removal methods:";
    for (std::size_t i = 0; i < n; ++i) {
        ss << " " << OutlierRemovalStrings[i];
        if(i != n - 1) ss << ",";
    }
    ss << ")";

    throw std::invalid_argument(ss.str());
}
