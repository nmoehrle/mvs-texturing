/*
 * Copyright (C) 2015, Nils Moehrle, Benjamin Richter
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <limits>

#include "icm_graph.h"

MRF_NAMESPACE_BEGIN

ICMGraph::ICMGraph(int num_sites, int) :
    sites(num_sites) {}

ENERGY_TYPE ICMGraph::compute_energy() {
    ENERGY_TYPE energy = 0;
    for (std::size_t i = 0; i < sites.size(); ++i) {
        Site const & site = sites[i];
        energy += site.data_cost + smooth_cost(i, site.label);
    }
    return energy;
}

ENERGY_TYPE ICMGraph::optimize(int num_iterations) {
    for (int i = 0; i < num_iterations; ++i) {
        for (std::size_t j = 0; j < sites.size(); ++j) {
            Site * site = &sites[j];
            /* Current cost */
            ENERGY_TYPE min_cost = std::numeric_limits<ENERGY_TYPE>::max(); //site->data_cost + smooth_cost(j, site->label);
            for (std::size_t k = 0; k < site->labels.size(); ++k) {
                ENERGY_TYPE cost = site->data_costs[k] + smooth_cost(j, site->labels[k]);
                if (cost < min_cost) {
                    min_cost = cost;
                    site->data_cost = site->data_costs[k];
                    site->label = site->labels[k];
                }
            }
        }
    }
    return compute_energy();
}

void ICMGraph::set_smooth_cost(SmoothCostFunction func) {
    smooth_cost_func = func;
}

void ICMGraph::set_neighbors(int site1, int site2) {
    sites[site1].neighbors.push_back(site2);
    sites[site2].neighbors.push_back(site1);
}


void ICMGraph::set_data_costs(int label, std::vector<SparseDataCost> const & costs) {
    for (std::size_t i = 0; i < costs.size(); ++i) {
        Site * site = &sites[costs[i].site];
        site->labels.push_back(label);
        int data_cost = costs[i].cost;
        site->data_costs.push_back(data_cost);

        if (data_cost < site->data_cost) {
            site->label = label;
            site->data_cost = data_cost;
        }
    }
}

int ICMGraph::what_label(int site) {
    return sites[site].label;
}

ENERGY_TYPE ICMGraph::smooth_cost(int site, int label) {
    ENERGY_TYPE smooth_cost = 0;
    for (int neighbor : sites[site].neighbors) {
         smooth_cost += smooth_cost_func(site, neighbor, label, sites[neighbor].label);
    }
    return smooth_cost;
}

int ICMGraph::num_sites() {
    return static_cast<int>(sites.size());
}

MRF_NAMESPACE_END
