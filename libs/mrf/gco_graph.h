/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef MRF_GCOGRAPH_HEADER
#define MRF_GCOGRAPH_HEADER

#include "graph.h"

#ifdef RESEARCH
#define GCO_ENERGYTYPE float
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

#endif /* RESEARCH */

#endif /* MRF_GCOGRAPH_HEADER */
