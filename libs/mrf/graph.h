/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef MRF_GRAPH_HEADER
#define MRF_GRAPH_HEADER

#include <vector>
#include <memory>

#define MRF_MAX_ENERGYTERM 10000000
#define MRF_NAMESPACE_BEGIN namespace mrf {
#define MRF_NAMESPACE_END }

MRF_NAMESPACE_BEGIN

typedef float ENERGY_TYPE;
typedef ENERGY_TYPE (*SmoothCostFunction)(int s1, int s2, int l1, int l2);

struct SparseDataCost {
    int site;
    ENERGY_TYPE cost;
};

enum SOLVER_TYPE {
    ICM,
    LBP,
    #ifdef RESEARCH
    GCO
    #endif
};

class Graph {
public:
    typedef std::shared_ptr<Graph> Ptr;
    virtual ~Graph() = default;
    virtual void set_smooth_cost(SmoothCostFunction func) = 0;
    virtual void set_data_costs(int label, std::vector<SparseDataCost> const & costs) = 0;
    virtual void set_neighbors(int site1, int site2) = 0;
    virtual ENERGY_TYPE compute_energy() = 0;
    virtual ENERGY_TYPE optimize(int num_iterations) = 0;
    virtual int what_label(int site) = 0;
    virtual int num_sites() = 0;

    static Graph::Ptr create(int num_sites, int num_labels, SOLVER_TYPE solver_type);
};

MRF_NAMESPACE_END

#endif /* MRF_GRAPH_HEADER */
