#include "MRF.h"
#include <algorithm>
#include <iostream>
#include <limits>

#ifdef RESEARCH
GCOWrapper::GCOWrapper(int num_sites, int num_lables) : gco(num_sites, num_lables) {
    /* GCoptimization uses rand() to create a random label order
     * - specify seed for repeadability. */
    srand(9313513);
    gco.setLabelOrder(true);
}

void GCOWrapper::set_smooth_cost(SmoothCostFunction func) {
    gco.setSmoothCost((GCoptimization::SmoothCostFn)func);
}


bool comp_spd_site(GCoptimization::SparseDataCost spd1, GCoptimization::SparseDataCost spd2){
    return spd1.site < spd2.site;
}

void GCOWrapper::set_data_costs(int label, std::vector<MRF::SparseDataCost> * costs) {
    costs = dynamic_cast<std::vector<GCoptimization::SparseDataCost>* >(costs);

    /* Sparse data costs must be sorted in increasing order of site ID */
    std::sort(costs->begin(), costs->end(), comp_spd_site);

    gco.setDataCost(label, costs->data(), costs->size());
}

MRF::ENERGY_TYPE GCOWrapper::compute_energy() {
    return static_cast<ENERGY_TYPE>(gco.compute_energy());
}

MRF::ENERGY_TYPE GCOWrapper::optimize(int num_iterations) {
    try {
        return static_cast<ENERGY_TYPE>(gco.expansion(num_iterations));
    } catch (GCException e) {
        std::cerr << e.message << std::endl;
        exit(EXIT_FAILURE);
    }
}

void GCOWrapper::set_neighbors(int site1, int site2) {
    gco.setNeighbors(site1, site2);
}

int GCOWrapper::what_label(int site) {
    return static_cast<int>(gco.whatLabel(site));
}
#endif

ICMSolver::ICMSolver(int num_sites, int) : sites(num_sites) {}

MRF::ENERGY_TYPE ICMSolver::compute_energy() {
    ENERGY_TYPE energy = 0;
    for (int i = 0; i < sites.size(); ++i) {
        Site const & site = sites[i];
        energy += site.data_cost + smooth_cost(i, site.label);
    }
    return energy;
}

MRF::ENERGY_TYPE ICMSolver::optimize(int num_iterations) {
    for (int i = 0; i < num_iterations; ++i) {
        for (int j = 0; j < sites.size(); ++j) {
            Site * site = &sites[j];
            /* Current cost */
            ENERGY_TYPE min_cost = std::numeric_limits<ENERGY_TYPE>::max(); //site->data_cost + smooth_cost(j, site->label);
            for (int k = 0; k < site->labels.size(); ++k) {
                ENERGY_TYPE cost = site->data_costs[k] + smooth_cost(j, site->labels[k]);
                if (cost < min_cost) {
                    min_cost = cost;
                    site->data_cost = site->data_costs[k];
                    site->label = site->labels[k];
                }
            }
        }
    }
    return compute_energy();
}

void ICMSolver::set_smooth_cost(MRF::SmoothCostFunction func) {
    smooth_cost_func = func;
}

void ICMSolver::set_neighbors(int site1, int site2) {
    sites[site1].neighbors.push_back(site2);
    sites[site2].neighbors.push_back(site1);
}


void ICMSolver::set_data_costs(int label, std::vector<MRF::SparseDataCost> * costs) {
    for (int i = 0; i < costs->size(); ++i) {
        Site * site = &sites[costs->at(i).site];
        site->labels.push_back(label);
        int data_cost = costs->at(i).cost;
        site->data_costs.push_back(data_cost);

        if (data_cost < site->data_cost) {
            site->label = label;
            site->data_cost = data_cost;
        }
    }
}

int ICMSolver::what_label(int site) {
    return sites[site].label;
}

MRF::ENERGY_TYPE ICMSolver::smooth_cost(int site, int label) {
    ENERGY_TYPE smooth_cost = 0;
    for (int neighbor : sites[site].neighbors) {
         smooth_cost += smooth_cost_func(site, neighbor, label, sites[neighbor].label);
    }
    return smooth_cost;
}

LBPSolver::LBPSolver(int num_sites, int) : vertices(num_sites) {}

MRF::ENERGY_TYPE LBPSolver::compute_energy() {
    ENERGY_TYPE energy = 0;

    #pragma omp parallel for reduction(+:energy)
    for (int vertex_idx = 0; vertex_idx < vertices.size(); ++vertex_idx) {
        Vertex const & vertex = vertices[vertex_idx];
        energy += vertex.data_cost;
    }

    #pragma omp parallel for reduction(+:energy)
    for (int edge_idx = 0; edge_idx < edges.size(); ++edge_idx) {
        DirectedEdge const & edge = edges[edge_idx];
        energy += smooth_cost_func(edge.v1, edge.v2, vertices[edge.v1].label, vertices[edge.v2].label);
    }

    return energy;
}

MRF::ENERGY_TYPE LBPSolver::optimize(int num_iterations) {
    for (int i = 0; i < num_iterations; ++i) {
        #pragma omp parallel for
        for (int edge_idx = 0; edge_idx < edges.size(); ++edge_idx) {
            DirectedEdge & edge = edges[edge_idx];
            std::vector<int> const & labels1 = vertices[edge.v1].labels;
            std::vector<int> const & labels2 = vertices[edge.v2].labels;
            for (int j = 0; j < labels2.size(); ++j) {
                int label2 = labels2[j];
                ENERGY_TYPE min_energy = std::numeric_limits<ENERGY_TYPE>::max();
                for (int k = 0; k < labels1.size(); ++k) {
                    int label1 = labels1[k];
                    ENERGY_TYPE energy = smooth_cost_func(edge.v1, edge.v2, label1, label2) + vertices[edge.v1].data_costs[k];

                    std::vector<int> const& incoming_edges1 = vertices[edge.v1].incoming_edges;
                    for (int n = 0; n < incoming_edges1.size(); ++n) {
                        DirectedEdge const& pre_edge = edges[incoming_edges1[n]];
                        if (pre_edge.v1 == edge.v2) continue;
                        energy += pre_edge.old_msg[k];
                    }

                    if (energy < min_energy)
                        min_energy = energy;
                }
                edge.new_msg[j] = min_energy;
            }
        }

        #pragma omp parallel for
        for (int edge_idx = 0; edge_idx < edges.size(); ++edge_idx) {
            DirectedEdge & edge = edges[edge_idx];
            edge.new_msg.swap(edge.old_msg);
            ENERGY_TYPE min_msg = std::numeric_limits<ENERGY_TYPE>::max();
            for (ENERGY_TYPE msg : edge.old_msg)
               min_msg = std::min(min_msg, msg);
            for (ENERGY_TYPE &msg : edge.old_msg)
               msg -= min_msg;
        }
    }

    #pragma omp parallel for
    for (int vertex_idx = 0; vertex_idx < vertices.size(); ++vertex_idx) {
        Vertex & vertex = vertices[vertex_idx];
        ENERGY_TYPE min_energy = std::numeric_limits<ENERGY_TYPE>::max();
        for (int j = 0; j < vertex.labels.size(); ++j) {
            ENERGY_TYPE energy = vertex.data_costs[j];
            for (int incoming_edge_idx : vertex.incoming_edges) {
                energy += edges[incoming_edge_idx].old_msg[j];
            }
            if (energy < min_energy) {
                min_energy = energy;
                vertex.label = vertex.labels[j];
                vertex.data_cost = vertex.data_costs[j];
            }
        }
    }

    return compute_energy();
}

void LBPSolver::set_smooth_cost(MRF::SmoothCostFunction func) {
    smooth_cost_func = func;
}

void LBPSolver::set_neighbors(int site1, int site2){
    edges.push_back(DirectedEdge(site1, site2));
    vertices[site2].incoming_edges.push_back(edges.size() - 1);
    edges.push_back(DirectedEdge(site2, site1));
    vertices[site1].incoming_edges.push_back(edges.size() - 1);
}


void LBPSolver::set_data_costs(int label, std::vector<MRF::SparseDataCost> * costs) {
    for (int i = 0; i < costs->size(); ++i) {
        Vertex * vertex = &vertices[costs->at(i).site];
        vertex->labels.push_back(label);
        int data_cost = costs->at(i).cost;
        vertex->data_costs.push_back(data_cost);

        if (data_cost < vertex->data_cost) {
            vertex->label = label;
            vertex->data_cost = data_cost;
        }

        for (int j : vertex->incoming_edges) {
            DirectedEdge &incoming_edge = edges[j];
            incoming_edge.old_msg.push_back(0);
            incoming_edge.new_msg.push_back(0);
        }
    }
}

int LBPSolver::what_label(int site) {
    return vertices[site].label;
}
