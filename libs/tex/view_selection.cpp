#include "view_selection.h"

void
fix_holes(UniGraph * graph) {
    /* Fix holes due to geometric inaccuracies. */
    std::vector<std::vector<std::size_t> > subgraphs;
    graph->get_subgraphs(0, &subgraphs);
    #pragma omp parallel for
    for (std::size_t i = 0; i < subgraphs.size(); ++i) {
        /* Only fix small holes, */
        const std::size_t max_hole_size = 10;
        if (subgraphs[i].size() < max_hole_size) {
            /* which are completely contained within another patch. */
            std::set<std::size_t> adj_labels;
            graph->get_adj_labels(subgraphs[i], &adj_labels);
            adj_labels.erase(0);
            if (adj_labels.size() == 1 && graph->get_min_conn(subgraphs[i]) == 3) {
                const std::size_t new_label = *(adj_labels.begin());
                for (std::size_t node : subgraphs[i]) {
                    graph->set_label(node, new_label);
                    //TODO check if projection is valid... (very likely)
                }
            }
        }
    }
}
