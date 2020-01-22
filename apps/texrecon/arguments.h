/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef ARGUMENTS_HEADER
#define ARGUMENTS_HEADER

#include "util/arguments.h"
#include "tex/settings.h"

/** Struct containing the commandline arguments at runtime. */
struct Arguments {
    std::string in_scene;
    std::string in_mesh;
    std::string out_prefix;

    std::string data_cost_file;
    std::string labeling_file;

    tex::Settings settings;

    bool write_timings;
    bool write_intermediate_results;
    bool write_view_selection_model;

    int num_threads;

    /** Returns a muliline string of the current arguments. */
    std::string to_string();
};

/**
 * Parses the commandline arguments.
 * @throws std::invalid_argument
 */
Arguments parse_args(int argc, char **argv);

#endif /* ARGUMENTS_HEADER */
