#pragma once

#include "Graph.h"

MRF_NAMESPACE_BEGIN

/** Implementation of the iterated conditional mode algrorithm. */
class ICMGraph : public Graph {
    private:
        struct Site {
            int label;
            int data_cost;
            std::vector<int> labels;
            std::vector<int> data_costs;
            std::vector<int> neighbors;
            Site() : label(0), data_cost(MRF_MAX_ENERGYTERM) {}
        };

        std::vector<Site> sites;
        SmoothCostFunction smooth_cost_func;
    public:
        ICMGraph(int num_sites, int num_labels);
        ENERGY_TYPE smooth_cost(int site, int label);

        void set_smooth_cost(SmoothCostFunction func);
        void set_data_costs(int label, std::vector<SparseDataCost> const & costs);
        void set_neighbors(int site1, int site2);
        ENERGY_TYPE compute_energy();
        ENERGY_TYPE optimize(int num_iterations);
        int what_label(int site);

        int num_sites();
};

MRF_NAMESPACE_END
