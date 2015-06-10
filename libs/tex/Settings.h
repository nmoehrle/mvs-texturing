#pragma once
#include "DataTerm.h"
#include "SmoothnessTerm.h"
#include "OutlierRemoval.h"

//TEX_NAMESPACE_BEGIN
struct Settings {
    bool verbose;

    DataTerm data_term;
    SmoothnessTerm smoothness_term;
    OutlierRemoval outlier_removal;

    bool geometric_visibility_test;
    bool global_seam_leveling;
    bool local_seam_leveling;
};
//TEX_NAMESPACE_END

