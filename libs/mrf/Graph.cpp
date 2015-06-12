#include "Graph.h"
#include "ICMGraph.h"
#include "LBPGraph.h"
#include "GCOGraph.h"

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
