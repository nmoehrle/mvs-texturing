/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include "icm_graph.h"
#include "lbp_graph.h"
#include "gco_graph.h"
#include "graph.h"

MRF_NAMESPACE_BEGIN

Graph::Ptr Graph::create(int num_sites, int num_labels, SOLVER_TYPE solver_type) {
    switch (solver_type) {
        case ICM: return Graph::Ptr(new ICMGraph(num_sites, num_labels));
        case LBP: return Graph::Ptr(new LBPGraph(num_sites, num_labels));
        #ifdef RESEARCH
        case GCO: return Graph::Ptr(new GCOGraph(num_sites, num_labels));
        #endif
        default: return Graph::Ptr(new LBPGraph(num_sites, num_labels));
    }
}

MRF_NAMESPACE_END
