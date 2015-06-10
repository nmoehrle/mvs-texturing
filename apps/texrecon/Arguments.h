#pragma once

#include "util/arguments.h"
#include "tex/Settings.h"

/** Struct containing the commandline arguments at runtime. */
struct Arguments {
    std::string in_scene;
    std::string in_mesh;
    std::string out_prefix;

    std::string data_cost_file;
    std::string labeling_file;

    Settings settings;

    bool write_timings;
    bool write_intermediate_results;
    bool write_view_selection_model;

    /** Returns a muliline string of the current arguments. */
    std::string to_string();
};

/**
 * Parses the commandline arguments.
 * @throws std::invalid_argument
 */
Arguments parse_args(int argc, char **argv);
