#pragma once
#include <vector>

#define MRF_MAX_ENERGYTERM 10000000

#ifdef RESEARCH
#include "GCoptimization.h"
#endif

class MRF {
    public:
        typedef long int ENERGY_TYPE;
        typedef int (*SmoothCostFunction)(int s1, int s2, int l1, int l2);

#ifdef RESEARCH
        typedef GCoptimization::SparseDataCost SparseDataCost;
#else
        struct SparseDataCost {
            int site;
            int cost;
        };
#endif

        virtual void set_smooth_cost(SmoothCostFunction func) = 0;
        virtual void set_data_costs(int label, std::vector<SparseDataCost>* costs) = 0;
        virtual void set_neighbors(int site1, int site2) = 0;
        virtual ENERGY_TYPE compute_energy() = 0;
        virtual ENERGY_TYPE optimize(int num_iterations) = 0;
        virtual int what_label(int site) = 0;
};

#ifdef RESEARCH
class GCOWrapper : public MRF {
    private:
        GCoptimizationGeneralGraph gco;

    public:
        GCOWrapper(int num_sites, int num_lables);

        void set_smooth_cost(SmoothCostFunction func);
        void set_data_costs(int label, std::vector<SparseDataCost>* costs);
        void set_neighbors(int site1, int site2);
        ENERGY_TYPE compute_energy();
        ENERGY_TYPE optimize(int num_iterations);
        int what_label(int site);
};
#endif

/** Implementation of the iterated conditional mode algrorithm. */
class ICMSolver : public MRF {
    private:
        struct Site {
            int label;
            int data_cost;
            std::vector<int> labels;
            std::vector<int> data_costs;
            std::vector<int> neighbors;
            Site() : label(0), data_cost(MRF_MAX_ENERGYTERM) {}
        };


        std::vector<Site> sites;
        SmoothCostFunction smooth_cost_func;
    public:
        ICMSolver(int num_sites, int num_labels);
        ENERGY_TYPE smooth_cost(int site, int label);

        void set_smooth_cost(SmoothCostFunction func);
        void set_data_costs(int label, std::vector<SparseDataCost>* costs);
        void set_neighbors(int site1, int site2);
        ENERGY_TYPE compute_energy();
        ENERGY_TYPE optimize(int num_iterations);
        int what_label(int site);
};

/** Implementation of the loopy belief propagation algorithm. */
class LBPSolver : public MRF {
    private:
        struct DirectedEdge {
            int v1;
            int v2;
            std::vector<ENERGY_TYPE> new_msg;
            std::vector<ENERGY_TYPE> old_msg;
            DirectedEdge(int v1, int v2) : v1(v1), v2(v2) {}
        };

        struct Vertex {
            int label;
            int data_cost;
            std::vector<int> labels;
            std::vector<int> data_costs;
            std::vector<int> incoming_edges;
            Vertex() : label(0), data_cost(MRF_MAX_ENERGYTERM) {}
        };

        std::vector<DirectedEdge> edges;
        std::vector<Vertex> vertices;
        SmoothCostFunction smooth_cost_func;
    public:
        LBPSolver(int num_sites, int num_labels);

        void set_smooth_cost(SmoothCostFunction func);
        void set_data_costs(int label, std::vector<SparseDataCost>* costs);
        void set_neighbors(int site1, int site2);
        ENERGY_TYPE compute_energy();
        ENERGY_TYPE optimize(int num_iterations);
        int what_label(int site);
};
