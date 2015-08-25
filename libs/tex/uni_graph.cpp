/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <limits>
#include <list>

#include "uni_graph.h"

UniGraph::UniGraph(std::size_t nodes) {
    adj_lists.resize(nodes);
    labels.resize(nodes);
    edges = 0;
}

void
UniGraph::get_subgraphs(std::size_t label,
    std::vector<std::vector<std::size_t> > * subgraphs) const {

    std::vector<bool> used(adj_lists.size(), false);

    for(std::size_t i = 0; i < adj_lists.size(); ++i) {
        if (labels[i] == label && !used[i]) {
            subgraphs->push_back(std::vector<std::size_t>());

            std::list<std::size_t> queue;

            queue.push_back(i);
            used[i] = true;

            while (!queue.empty()) {
                std::size_t node = queue.front();
                queue.pop_front();

                subgraphs->back().push_back(node);

                /* Add all unused neighbours with the same label to the queue. */
                std::vector<std::size_t> const & adj_list = adj_lists[node];
                for(std::size_t j = 0; j < adj_list.size(); ++j) {
                    std::size_t adj_node = adj_list[j];
                    assert(adj_node < labels.size() && adj_node < used.size());
                    if (labels[adj_node] == label && !used[adj_node]){
                        queue.push_back(adj_node);
                        used[adj_node] = true;
                    }
                }
            }
        }
    }
}
