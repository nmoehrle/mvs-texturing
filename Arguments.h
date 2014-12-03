#pragma once

#include "util/arguments.h"
#include "DataTerm.h"
#include "SmoothnessTerm.h"
#include "OutlierRemoval.h"

/** Struct containing the commandline arguments at runtime. */
struct Arguments {
    std::string in_scene;
    std::string in_mesh;
    std::string out_prefix;

    std::string data_cost_file;
    std::string labeling_file;

    DataTerm data_term;
    SmoothnessTerm smoothness_term;
    OutlierRemoval outlier_removal;

    bool write_view_selection_model;
    bool write_data_term_histograms;
    bool write_mrf_energies;
    bool global_seam_leveling;
    bool local_seam_leveling;

    /** Returns a muliline string of the current arguments. */
    std::string to_string();
};

/**
 * Parses the commandline arguments.
 * @throws std::invalid_argument
 */
Arguments parse_args(int argc, char **argv);

/**
 * Parses the string s to its corresponding SmoothnessTerm.
 * @throws std::invalid_argument
 */
SmoothnessTerm parse_smoothness_term(std::string s);

/**
 * Parses the string s to its corresponding DataTerm.
 * @throws std::invalid_argument
 */
DataTerm parse_data_term(std::string s);

/**
 * @param s
 * @return The outlier removal method corresponding to s.
 */
OutlierRemoval parse_outlier_removal(std::string s);
