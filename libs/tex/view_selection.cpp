#include "texturing.h"
#include "util.h"

TEX_NAMESPACE_BEGIN

bool IGNORE_LUMINANCE = false;

/** Potts model */
float
potts(int, int, int l1, int l2) {
    return l1 == l2 && l1 != 0 && l2 != 0 ? 0 : 1 * MRF_MAX_ENERGYTERM;
}

struct FaceInfo {
    std::size_t component;
    std::size_t id;
};

/** Setup the neighborhood of the MRF. */
void
set_neighbors(UniGraph const & graph, std::vector<FaceInfo> const & face_infos,
    std::vector<mrf::Graph::Ptr> const & mrfs) {
    for (std::size_t i = 0; i < graph.num_nodes(); ++i) {
        std::vector<std::size_t> adj_faces = graph.get_adj_nodes(i);
        for (std::size_t j = 0; j < adj_faces.size(); ++j) {
            std::size_t adj_face = adj_faces[j];
            /* The solver expects only one call of setNeighbours for two neighbours a and b. */
            if (i < adj_face) {
                assert(face_infos[i].component == face_infos[adj_face].component);
                const std::size_t component = face_infos[i].component;
                const std::size_t cid1 = face_infos[i].id;
                const std::size_t cid2 = face_infos[adj_face].id;
                mrfs[component]->set_neighbors(cid1, cid2);
            }
        }
    }
}

/** Set the data costs of the MRF. */
void
set_data_costs(std::vector<FaceInfo> const & face_infos, ST const & data_costs,
    std::vector<mrf::Graph::Ptr> const & mrfs) {

    /* Set data costs for all labels except label 0 (undefined) */
    for (std::size_t i = 0; i < data_costs.rows(); i++) {
        ST::Row const & data_costs_for_label = data_costs.row(i);

        std::vector<std::vector<mrf::SparseDataCost> > costs(mrfs.size());
        for(std::size_t j = 0; j < data_costs_for_label.size(); j++) {
            const std::size_t id = data_costs_for_label[j].first;
            const float data_cost = data_costs_for_label[j].second;
            const std::size_t component = face_infos[id].component;
            const std::size_t cid = face_infos[id].id;
            //TODO change index type of mrf::Graph
            costs[component].push_back({static_cast<int>(cid), data_cost});
        }

        int label = i + 1;

        for (std::size_t j = 0; j < mrfs.size(); ++j) {
            mrfs[j]->set_data_costs(label, costs[j]);
        }
    }

    for (std::size_t i = 0; i < mrfs.size(); ++i) {
        /* Set costs for undefined label */
        std::vector<mrf::SparseDataCost> costs(mrfs[i]->num_sites());
        for (std::size_t j = 0; j < costs.size(); j++) {
            costs[j].site = j;
            costs[j].cost = MRF_MAX_ENERGYTERM;
        }
        mrfs[i]->set_data_costs(0, costs);
    }
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
optimize(mrf::Graph * mrf, std::size_t component, bool verbose = false) {
    util::WallTimer timer;

    std::vector<mrf::ENERGY_TYPE> energies;

    mrf::ENERGY_TYPE const zero = mrf::ENERGY_TYPE(0);
    mrf::ENERGY_TYPE last_energy = zero;
    mrf::ENERGY_TYPE energy = mrf->compute_energy();
    mrf::ENERGY_TYPE diff = last_energy - energy;
    unsigned int i = 0;

    #ifdef RESEARCH
    #ifdef _OPENMP
    if (component == 0) {
        std::cout << "\t" << "Parallel optimization with " <<
            omp_get_num_threads() << " threads" << std::endl;
        std::cout << "\tComp\tIter\tEnergy\t\tRuntime" << std::endl;
    }
    #else
    if (verbose) std::cout << "\tComp\tIter\tEnergy\t\tRuntime" << std::endl;
    #endif
    #else
    if (verbose) std::cout << "\tComp\tIter\tEnergy\t\tRuntime" << std::endl;
    #endif
    while (diff != zero) {
        #pragma omp critical
        if (verbose) {
            std::cout << "\t" << component << "\t" << i << "\t" << energy
                << "\t" << timer.get_elapsed_sec() << std::endl;
        }
        energies.push_back(energy);
        last_energy = energy;
        i++;
        energy = mrf->optimize(1);
        diff = last_energy - energy;
        if (diff <= 0) break;
    }

    #pragma omp critical
    if (verbose) {
        std::cout << "\t" << component << "\t" << i << "\t" << energy << std::endl;
        if (diff == zero) {
            std::cout << "\t" << component << "\t" << "Converged" << std::endl;
        }
        if (diff < zero) {
            std::cout << "\t" << component << "\t"
                << "Increase of energy detected! Aborting..." << std::endl;
        }
    }
    //if (conf.write_mrf_energies)
    //    write_vector_to_csv(conf.out_prefix + "_mrf_energies.csv", energies, "Energy");

}

void
view_selection(ST const & data_costs, UniGraph * graph, Settings const & settings) {
    UniGraph mgraph(*graph);
    isolate_unseen_faces(&mgraph, data_costs);

    std::vector<FaceInfo> face_infos(mgraph.num_nodes());
    std::vector<std::vector<std::size_t> > components;
    mgraph.get_subgraphs(0, &components);
    for (std::size_t i = 0; i < components.size(); ++i) {
        for (std::size_t j = 0; j < components[i].size(); ++j) {
            face_infos[components[i][j]] = {i, j};
        }
    }

    #ifdef RESEARCH
    mrf::SOLVER_TYPE solver_type = mrf::GCO;
    #else
    mrf::SOLVER_TYPE solver_type = mrf::LBP;
    #endif

    /* Label 0 is undefined. */
    const std::size_t num_labels = data_costs.rows() + 1;
    std::vector<mrf::Graph::Ptr> mrfs(components.size());
    for (std::size_t i = 0; i < components.size(); ++i) {
        mrfs[i] = mrf::Graph::create(components[i].size(), num_labels, solver_type);
    }

    /* Set neighbors must be called prior to set_data_costs (LBP). */
    set_neighbors(mgraph, face_infos, mrfs);

    set_data_costs(face_infos, data_costs, mrfs);

    #ifdef RESEARCH
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (std::size_t i = 0; i < components.size(); ++i) {
        switch (settings.smoothness_term) {
            case POTTS:
                mrfs[i]->set_smooth_cost(*potts);
            break;
        }

        optimize(mrfs[i].get(), i, mrfs[i]->num_sites() > 10000);

        /* Extract resulting labeling from MRF. */
        for (std::size_t j = 0; j < components[i].size(); ++j) {
            int label = mrfs[i]->what_label(static_cast<int>(j));
            assert(0 <= label && static_cast<std::size_t>(label) < num_labels);
            graph->set_label(components[i][j], static_cast<std::size_t>(label));
        }
    }
}

TEX_NAMESPACE_END
