#include "texturing.h"
#include "util.h"

bool IGNORE_LUMINANCE = false;

/** Potts model */
float
potts(int s1, int s2, int l1, int l2) {
    /* Suppress compiler warning because of unused variable. */
    (void) s1; (void) s2;
    return l1 == l2 && l1 != 0 && l2 != 0 ? 0 : 1 * MRF_MAX_ENERGYTERM;
}

/** Struct representing a non directed edge. */
struct Edge {
    std::size_t v1_id;
    std::size_t v2_id;

    Edge(std::size_t _v1_id, std::size_t _v2_id) {
        if(_v1_id < _v2_id) {
            v1_id = _v1_id;
            v2_id = _v2_id;
        } else {
            v1_id = _v2_id;
            v2_id = _v1_id;
        }
    }
};

/** Lexical comparison of non directed edges. */
bool
operator<(Edge const & edge1, Edge const & edge2) {
    return edge1.v1_id < edge2.v1_id || (edge1.v1_id == edge2.v1_id && edge1.v2_id < edge2.v2_id);
}

/** Setup the neighborhood of the MRF. */
void
set_neighbors(mrf::Graph * mrf, UniGraph const & graph) {
    for (std::size_t i = 0; i < graph.num_nodes(); ++i) {
        std::vector<std::size_t> adj_faces = graph.get_adj_nodes(i);
        for (std::size_t j = 0; j < adj_faces.size(); ++j) {
            std::size_t adj_face = adj_faces[j];
            /* The solver expects only one call of setNeighbours for two neighbours a and b. */
            if (i < adj_face) mrf->set_neighbors(i, adj_face);
        }
    }
}



/** Set the data costs of the MRF. */
void
set_data_costs(mrf::Graph * mrf, ST const & data_costs){
    /* Set data costs for all labels except label 0 (undefined) */
    for (std::size_t i = 0; i < data_costs.rows(); i++) {
        ST::Row const & data_costs_for_label = data_costs.row(i);

        std::vector<mrf::SparseDataCost> costs(data_costs_for_label.size());
        for(std::size_t j = 0; j < costs.size(); j++) {
            costs[j].site = data_costs_for_label[j].first;
            costs[j].cost = data_costs_for_label[j].second;
        }


        int label = i + 1;
        mrf->set_data_costs(label, costs);
    }

    /* Set costs for undefined label */
    std::vector<mrf::SparseDataCost> costs(data_costs.cols());
    for (std::size_t i = 0; i < costs.size(); i++) {
        costs[i].site = i;
        costs[i].cost = MRF_MAX_ENERGYTERM;
    }
    mrf->set_data_costs(0, costs);
}

/** Remove all edges of nodes which corresponding face has not been seen in any texture view. */
void
isolate_unseen_faces(UniGraph * graph, ST const & data_costs) {
    int num_unseen_faces = 0;
    for (std::uint32_t i = 0; i < data_costs.cols(); i++) {
        ST::Column const & data_costs_for_face = data_costs.col(i);

        if (data_costs_for_face.size() == 0) {
            num_unseen_faces++;

            std::vector<std::size_t> const & adj_nodes = graph->get_adj_nodes(i);
            for (std::size_t j = 0; j < adj_nodes.size(); j++)
                graph->remove_edge(i, adj_nodes[j]);
        }

    }
    std::cout << "\t" << num_unseen_faces << " faces have not been seen by a view." << std::endl;
}

void
build_mrf(UniGraph graph, ST const & data_costs, mrf::Graph * mrf, Settings const & settings) {
    /* Graph is copied...*/
    isolate_unseen_faces(&graph, data_costs);

    /* Set neighbors must be called prior to set_data_costs (LBP). */
    set_neighbors(mrf, graph);

    set_data_costs(mrf, data_costs);

    switch (settings.smoothness_term) {
        case EDI:
            std::cerr << "\tThe smoothness term EDI has been removed - falling back to POTTS" << std::endl;
        case POTTS:
            mrf->set_smooth_cost(*potts);
        break;
    }
}
