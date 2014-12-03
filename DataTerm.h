#pragma once

#include <sstream>

/** Enum representing a data term. */
enum DataTerm {
    AREA = 0,
    GMI = 1
};

char const * const DataTermStrings [] = {"area", "gmi"};

inline
DataTerm parse_data_term(std::string s){
    std::size_t n = sizeof(DataTermStrings) / sizeof(char const *);
    for (std::size_t i = 0; i < n; ++i) {
        if (s == DataTermStrings[i]){
            return static_cast<DataTerm>(i);
        }
    }

    std::stringstream ss;
    ss << "Invalid data term: " << s << " (Available data terms:";
    for (std::size_t i = 0; i < n; ++i) {
        ss << " " << DataTermStrings[i];
        if(i != n - 1) ss << ",";
    }
    ss << ")";

    throw std::invalid_argument(ss.str());
}
