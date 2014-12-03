#pragma once

#include <sstream>

/** Enum representing a smoothness term. */
enum SmoothnessTerm {
    POTTS = 0,
    EDI = 1
};

char const * const SmoothnessTermStrings [] = {"potts", "edi"};

inline
SmoothnessTerm parse_smoothness_term(std::string s) {
    std::size_t n = sizeof(SmoothnessTermStrings) / sizeof(char const *);
    for (std::size_t i = 0; i < n; ++i) {
        if (s == SmoothnessTermStrings[i]) {
            return static_cast<SmoothnessTerm>(i);
        }
    }


    std::stringstream ss;
    ss << "Invalid smoothness term: " << s << " (Available smoothness terms:";
    for (std::size_t i = 0; i < n; ++i) {
        ss << " " << SmoothnessTermStrings[i];
        if(i != n - 1) ss << ",";
    }
    ss << ")";

    throw std::invalid_argument(ss.str());
}
