#pragma once

#include "MRF.h"

/** Implementation of the loopy belief propagation algorithm. */
class LBPSolver : public MRF {
    private:
        struct DirectedEdge {
            int v1;
            int v2;
            std::vector<ENERGY_TYPE> new_msg;
            std::vector<ENERGY_TYPE> old_msg;
            DirectedEdge(int v1, int v2) : v1(v1), v2(v2) {}
        };

        struct Vertex {
            int label;
            int data_cost;
            std::vector<int> labels;
            std::vector<int> data_costs;
            std::vector<int> incoming_edges;
            Vertex() : label(0), data_cost(MRF_MAX_ENERGYTERM) {}
        };

        std::vector<DirectedEdge> edges;
        std::vector<Vertex> vertices;
        SmoothCostFunction smooth_cost_func;
    public:
        LBPSolver(int num_sites, int num_labels);

        void set_smooth_cost(SmoothCostFunction func);
        void set_data_costs(int label, std::vector<SparseDataCost> const & costs);
        void set_neighbors(int site1, int site2);
        ENERGY_TYPE compute_energy();
        ENERGY_TYPE optimize(int num_iterations);
        int what_label(int site);
};
