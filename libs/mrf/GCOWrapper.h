#pragma once

#include "MRF.h"

#ifdef RESEARCH
#define GCO_ENERGYTYPE float
#include "GCoptimization.h"

class GCOWrapper : public MRF {
    private:
        GCoptimizationGeneralGraph gco;

    public:
        GCOWrapper(int num_sites, int num_lables);

        void set_smooth_cost(SmoothCostFunction func);
        void set_data_costs(int label, std::vector<SparseDataCost> const & costs);
        void set_neighbors(int site1, int site2);
        ENERGY_TYPE compute_energy();
        ENERGY_TYPE optimize(int num_iterations);
        int what_label(int site);
};
#endif
