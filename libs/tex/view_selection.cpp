/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <util/timer.h>

#include "util.h"
#include "texturing.h"
#include "mapmap/full.h"

TEX_NAMESPACE_BEGIN

void
view_selection(DataCosts const & data_costs, UniGraph * graph, Settings const &) {
    using uint_t = unsigned int;
    using cost_t = float;
    constexpr uint_t simd_w = mapmap::sys_max_simd_width<cost_t>();
    using unary_t = mapmap::UnaryTable<cost_t, simd_w>;
    using pairwise_t = mapmap::PairwisePotts<cost_t, simd_w>;

    /* Construct graph */
    mapmap::Graph<cost_t> mgraph(graph->num_nodes());

    for (std::size_t i = 0; i < graph->num_nodes(); ++i) {
        if (data_costs.col(i).empty()) continue;

        std::vector<std::size_t> adj_faces = graph->get_adj_nodes(i);
        for (std::size_t j = 0; j < adj_faces.size(); ++j) {
            std::size_t adj_face = adj_faces[j];
            if (data_costs.col(adj_face).empty()) continue;

            /* Uni directional */
            if (i < adj_face) {
                mgraph.add_edge(i, adj_face, 1.0f);
            }
        }
    }
    mgraph.update_components();

    mapmap::LabelSet<cost_t, simd_w> label_set(graph->num_nodes(), false);
    for (std::size_t i = 0; i < data_costs.cols(); ++i) {
        DataCosts::Column const & data_costs_for_node = data_costs.col(i);

        std::vector<mapmap::_iv_st<cost_t, simd_w> > labels;
        if (data_costs_for_node.empty()) {
            labels.push_back(0);
        } else {
            labels.resize(data_costs_for_node.size());
            for(std::size_t j = 0; j < data_costs_for_node.size(); ++j) {
                labels[j] = data_costs_for_node[j].first + 1;
            }
        }

        label_set.set_label_set_for_node(i, labels);
    }

    std::vector<unary_t> unaries;
    unaries.reserve(data_costs.cols());
    pairwise_t pairwise(1.0f);
    for (std::size_t i = 0; i < data_costs.cols(); ++i) {
        DataCosts::Column const & data_costs_for_node = data_costs.col(i);

        std::vector<mapmap::_s_t<cost_t, simd_w> > costs;
        if (data_costs_for_node.empty()) {
            costs.push_back(1.0f);
        } else {
            costs.resize(data_costs_for_node.size());
            for(std::size_t j = 0; j < data_costs_for_node.size(); ++j) {
                float cost = data_costs_for_node[j].second;
                costs[j] = cost;
            }

        }

        unaries.emplace_back(i, &label_set);
        unaries.back().set_costs(costs);
    }

    mapmap::StopWhenReturnsDiminish<cost_t, simd_w> terminate(5, 0.01);
    std::vector<mapmap::_iv_st<cost_t, simd_w> > solution;

    auto display = [](const mapmap::luint_t time_ms,
            const mapmap::_iv_st<cost_t, simd_w> objective) {
        std::cout << "\t\t" << time_ms / 1000 << "\t" << objective << std::endl;
    };

    /* Create mapMAP solver object. */
    mapmap::mapMAP<cost_t, simd_w> solver;
    solver.set_graph(&mgraph);
    solver.set_label_set(&label_set);
    for(std::size_t i = 0; i < graph->num_nodes(); ++i)
        solver.set_unary(i, &unaries[i]);
    solver.set_pairwise(&pairwise);
    solver.set_logging_callback(display);
    solver.set_termination_criterion(&terminate);

    /* Pass configuration arguments (optional) for solve. */
    mapmap::mapMAP_control ctr;
    ctr.use_multilevel = true;
    ctr.use_spanning_tree = true;
    ctr.use_acyclic = true;
    ctr.spanning_tree_multilevel_after_n_iterations = 5;
    ctr.force_acyclic = true;
    ctr.min_acyclic_iterations = 5;
    ctr.relax_acyclic_maximal = true;
    ctr.tree_algorithm = mapmap::LOCK_FREE_TREE_SAMPLER;

    /* Set true for deterministic (but slower) mapMAP execution. */
    ctr.sample_deterministic = false;
    ctr.initial_seed = 548923723;

    std::cout << "\tOptimizing:\n\t\tTime[s]\tEnergy" << std::endl;
    solver.optimize(solution, ctr);

    /* Label 0 is undefined. */
    std::size_t num_labels = data_costs.rows() + 1;
    std::size_t undefined = 0;
    /* Extract resulting labeling from solver. */
    for (std::size_t i = 0; i < graph->num_nodes(); ++i) {
        int label = label_set.label_from_offset(i, solution[i]);
        if (label < 0 || num_labels <= static_cast<std::size_t>(label)) {
            throw std::runtime_error("Incorrect labeling");
        }
        if (label == 0) undefined += 1;
        graph->set_label(i, static_cast<std::size_t>(label));
    }
    std::cout << '\t' << undefined << " faces have not been seen" << std::endl;
}

TEX_NAMESPACE_END
