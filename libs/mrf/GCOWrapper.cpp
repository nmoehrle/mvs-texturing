#include "GCOWrapper.h"

#include <iostream>
#include <algorithm>

#ifdef RESEARCH
GCOWrapper::GCOWrapper(int num_sites, int num_lables) : gco(num_sites, num_lables) {
    /* GCoptimization uses rand() to create a random label order
     * - specify seed for repeadability. */
    srand(9313513);
    gco.setLabelOrder(true);
}

void GCOWrapper::set_smooth_cost(SmoothCostFunction func) {
    gco.setSmoothCost((GCoptimization::SmoothCostFn) func);
}


bool comp_spd_site(GCoptimization::SparseDataCost spd1, GCoptimization::SparseDataCost spd2) {
    return spd1.site < spd2.site;
}

void GCOWrapper::set_data_costs(int label, std::vector<SparseDataCost> const & costs) {
    /* Sparse data costs must be sorted in increasing order of site ID */
    std::vector<GCoptimization::SparseDataCost> gcocosts(costs.size());
    for (std::size_t i = 0; i < costs.size(); ++i) {
        gcocosts[i] = {costs[i].site, costs[i].cost};
    }
    std::sort(gcocosts.begin(), gcocosts.end(), comp_spd_site);
    gco.setDataCost(label, gcocosts.data(), gcocosts.size());
}

ENERGY_TYPE GCOWrapper::compute_energy() {
    return static_cast<ENERGY_TYPE>(gco.compute_energy());
}

ENERGY_TYPE GCOWrapper::optimize(int num_iterations) {
    try {
        return static_cast<ENERGY_TYPE>(gco.expansion(num_iterations));
    } catch (GCException e) {
        std::cerr << e.message << std::endl;
        exit(EXIT_FAILURE);
    }
}

void GCOWrapper::set_neighbors(int site1, int site2) {
    gco.setNeighbors(site1, site2);
}

int GCOWrapper::what_label(int site) {
    return static_cast<int>(gco.whatLabel(site));
}
#endif

//MRF_NAMESPACE_END
