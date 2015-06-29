#pragma once

#include "Graph.h"

#ifdef RESEARCH
#include "GCoptimization.h"

MRF_NAMESPACE_BEGIN

class GCOGraph : public Graph {
    private:
        GCoptimizationGeneralGraph gco;

    public:
        GCOGraph(int num_sites, int num_lables);

        void set_smooth_cost(SmoothCostFunction func);
        void set_data_costs(int label, std::vector<SparseDataCost> const & costs);
        void set_neighbors(int site1, int site2);
        ENERGY_TYPE compute_energy();
        ENERGY_TYPE optimize(int num_iterations);
        int what_label(int site);

        int num_sites();
        int num_labels();
};

MRF_NAMESPACE_END

#endif
