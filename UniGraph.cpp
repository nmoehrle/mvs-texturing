#include "UniGraph.h"
#include <limits>

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

void
UniGraph::get_adj_labels(std::vector<std::size_t> const & nodes,
    std::set<std::size_t> * adj_labels) const {

    for (std::size_t node : nodes) {
        std::vector<std::size_t> adj_nodes = get_adj_nodes(node);
        for (std::size_t adj_node : adj_nodes) {
            adj_labels->insert(get_label(adj_node));
        }
    }
}

std::size_t UniGraph::get_min_conn(std::vector<std::size_t> const & nodes) const {
    std::size_t min_conn = std::numeric_limits<std::size_t>::max();
    for (std::size_t node : nodes) {
        min_conn = std::min(min_conn, get_adj_nodes(node).size());
    }
    return min_conn;
}
