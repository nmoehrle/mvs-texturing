#include "MRF.h"

#include "ICMSolver.h"
#include "LBPSolver.h"
#include "GCOWrapper.h"
MRF::Ptr MRF::create(int num_sites, int num_labels, SOLVER_TYPE solver_type) {
    switch (solver_type) {
        case ICM: return MRF::Ptr(new ICMSolver(num_sites, num_labels));
        case LBP: return MRF::Ptr(new LBPSolver(num_sites, num_labels));
        #ifdef RESEARCH
        case GCO: return MRF::Ptr(new GCOWrapper(num_sites, num_labels));
        #endif
    }
}
