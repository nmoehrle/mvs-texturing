/*
 * Copyright (C) 2015, Michael Waechter
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <vector>
#include <map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <random>
#include <iterator>

#include "uni_graph.h"

/**
 * Finds connected sets of vertices within the given graph and gives all connected vertices the same label.
 *
 * @param graph
 * @return the number of connected components
 */
std::size_t
label_connected_subgraphs(UniGraph * graph)
{
    std::size_t const num_nodes = graph->num_nodes();

    /* Set all nodes to an invalid label so we can detect this later. */
    for (std::size_t node = 0; node < num_nodes; ++node)
        graph->set_label(node, num_nodes + 1);

    std::vector<bool> visited(num_nodes, false);

    std::size_t num_components = 0;

    auto first_unvisited_node = visited.begin();

    /* While not all nodes have been visited yet. */
    while ((first_unvisited_node = std::find(first_unvisited_node, visited.end(), false)) != visited.end()) {

        /* Put first unvisited node into queue to start a new flooding algorithm from there. */
        std::vector<std::size_t> queue;
        queue.push_back(*first_unvisited_node);

        /* Flooding algorithm */
        while (!queue.empty()) {
            std::size_t const current_node = queue.back();
            queue.pop_back();
            visited.at(current_node) = true;
            graph->set_label(current_node, num_components);
            std::vector<std::size_t> const & neighbors = graph->get_adj_nodes(current_node);
            for (std::size_t const neighbor : neighbors)
                if (!visited.at(neighbor))
                    queue.push_back(neighbor);
        }

        /* after the queue is empty all nodes in this connected component
         * have been found and we increase the component count */
        ++num_components;
    }

    for (std::size_t node = 0; node < num_nodes; ++node)
        assert(graph->get_label(node) < num_components);

    return num_components;
}

/**
 * @param graph
 * @param num_labels the number of different labels in this graph
 * @return the number of occurences of the label that occurs most often in the graph
 */
std::size_t occurences_of_most_frequent_label(UniGraph const * graph) {
    std::map<std::size_t, std::size_t> label_frequencies;
    for (std::size_t node = 0; node < graph->num_nodes(); ++node)
        label_frequencies[graph->get_label(node)]++;

    const std::size_t most_frequent_label = std::max_element(label_frequencies.begin(), label_frequencies.end(), label_frequencies.value_comp())->second;
    return label_frequencies.at(most_frequent_label);
}

/**
 * Takes a subgraph of the given graph (all nodes in the graph with the given label),
 * partitions this subgraph into even smaller subgraphs (using something similar to k-means),
 * and gives all small subgraphs a unique label (using the given min_label).
 *
 * @param graph
 * @param label_of_connected_component
 * @param size_of_largest_partition
 * @param min_label_for_partition_labeling
 * @return the number of generated partitions
 */
std::size_t
partition_connected_component(UniGraph * graph, std::size_t label_of_connected_component, std::size_t partition_size, std::size_t min_label_for_partition_labeling)
{
    typedef std::size_t Node;
    typedef std::size_t Label;
    typedef std::vector<Node> Nodes;

    Nodes nodes;
    for (Node node = 0; node < graph->num_nodes(); ++node)
        if (graph->get_label(node) == label_of_connected_component)
            nodes.push_back(node);

    const std::size_t num_partitions = (nodes.size() + partition_size - 1) / partition_size; // division and rounding up

    /********* k-means clustering *******/

    const std::size_t num_kmeans_iterations = 100;

    Nodes centroids;

    /* Draw centroids randomly. */
    std::default_random_engine generator;
    std::uniform_int_distribution<std::size_t> distribution(0, nodes.size() - 1);
    for(std::size_t partition = 0; partition < num_partitions; ++partition) {
        Node centroid = std::numeric_limits<Node>::max();
        while (std::find(centroids.begin(), centroids.end(), centroid) != centroids.end())
            centroid = nodes.at(distribution(generator));
        centroids.push_back(centroid);
    }

    for (std::size_t kmeans_iteration = 0; kmeans_iteration < num_kmeans_iterations; ++kmeans_iteration) {
        const Label unvisited = std::numeric_limits<Label>::max();
        for (Node const & node : nodes)
            graph->set_label(node, unvisited);

        /* Put centroids into queues. */
        std::vector<Nodes> queues(num_partitions);
        for (std::size_t i = 0; i < num_partitions; ++i)
            queues.at(i).push_back(centroids.at(i));

        /* Grow regions starting from centroids */
        while (std::any_of(queues.begin(), queues.end(), [](Nodes const & queue){return !queue.empty();})) {
#pragma omp parallel for
            for (std::size_t queue_id = 0; queue_id < queues.size(); ++queue_id) {
                Nodes & old_queue = queues.at(queue_id);
                std::unordered_set<Node> new_queue;
                for (Node node : old_queue)
                    graph->set_label(node, min_label_for_partition_labeling + queue_id); // there is a race condition for partition boundary nodes but we don't care
                for (Node node : old_queue) {
                    /* Copy all unvisited (and not yet inserted) neighbors into new queue. */
                    for (Node neighbor : graph->get_adj_nodes(node))
                        if (graph->get_label(neighbor) == unvisited)
                            new_queue.insert(neighbor);
                }

                old_queue.clear();
                old_queue.insert(old_queue.begin(), new_queue.begin(), new_queue.end());
            }
        }

        /* If we are in the final iteration we stop here to keep the graph labels
         * (they would be removed in the following region shrinking step). */
        if (kmeans_iteration == num_kmeans_iterations - 1)
            break;

        /* Put partition boundary nodes into queues. */
        for (Node const node : nodes) {
            Label const cur_label = graph->get_label(node);
            std::size_t const cur_queue = cur_label - min_label_for_partition_labeling;
            Nodes const & neighbors = graph->get_adj_nodes(node);
            /* Each node, where any of its neighbors has a different label, is a boundary node. */
            if (std::any_of(neighbors.begin(), neighbors.end(), [graph, cur_label]
                (Node const neighbor) { return graph->get_label(neighbor) != cur_label; } ))
                queues.at(cur_queue).push_back(node);
        }

        /* Shrink regions starting from boundaries to obtain new centroids. */
#pragma omp parallel for
        for (std::size_t queue_id = 0; queue_id < queues.size(); ++queue_id) {
            Nodes & old_queue = queues.at(queue_id);
            while (!old_queue.empty()){
                std::unordered_set<Node> new_queue;
                for (Node node : old_queue)
                    graph->set_label(node, unvisited);
                for (Node node : old_queue) {
                    /* Copy all neighbors that have not yet been marked (and have not yet been inserted) into new queue. */
                    for (Node neighbor : graph->get_adj_nodes(node))
                        if (graph->get_label(neighbor) == min_label_for_partition_labeling + queue_id)
                            new_queue.insert(neighbor);
                }

                /* If the new queue is empty we are (almost) finished and use a random node from the old queue as new centroid. */
                if (new_queue.empty()) {
                    std::uniform_int_distribution<std::size_t> distribution(0, old_queue.size() - 1);
                    centroids.at(queue_id) = old_queue.at(distribution(generator));
                }

                /* Replace old queue with new one. */
                old_queue.clear();
                old_queue.insert(old_queue.begin(), new_queue.begin(), new_queue.end());
            }
        }
    }

    return num_partitions;
}

/**
 * This is best explained by an example:
 *
 * Given the following graph: 8 nodes, 3 nodes have label 0, 4 nodes have label 1, and 1 node has label 2.<br>
 * First the labels are sorted by relative frequency: label 1 = 4 occurences, label 0 = 3 occurences, label 2 = 1 occurence<br>
 * Second the labels are reassigned: label 1 --> 0, label 0 --> 1, label 2 --> 2.
 *
 * @param graph
 */
void
reassign_labels_by_relative_frequency(UniGraph * graph)
{
    /* Compute frequencies. */
    std::map<std::size_t, std::size_t> label_frequencies; // map with key = label and value = frequency
    for (std::size_t node = 0; node < graph->num_nodes(); ++node)
        label_frequencies[graph->get_label(node)]++;

    /* Insert elements into vector in a sorted way. */
    std::vector<std::size_t> sorted_labels;
    while (!label_frequencies.empty()) {
        auto iterator_to_largest_label = std::max_element(label_frequencies.begin(), label_frequencies.end(), label_frequencies.value_comp());
        sorted_labels.push_back(iterator_to_largest_label->first);
        label_frequencies.erase(iterator_to_largest_label);
    }

    /* Build map from old to new labels. */
    std::map<std::size_t, std::size_t> label_mapping;
    for (std::size_t new_label = 0; new_label < sorted_labels.size(); ++new_label)
        label_mapping[sorted_labels.at(new_label)] = new_label;

    /* Relabel all nodes. */
    for (std::size_t node = 0; node < graph->num_nodes(); ++node)
        graph->set_label(node, label_mapping.at(graph->get_label(node)));
}

void
partition_mesh(UniGraph *graph, std::size_t min_num_partitions)
{
    const std::size_t num_submeshes = label_connected_subgraphs(graph);
    const std::size_t largest_submesh_size = occurences_of_most_frequent_label(graph);
    const std::size_t partition_size = (largest_submesh_size + min_num_partitions - 1) / min_num_partitions; // division and rounding up

    std::size_t min_unused_label = num_submeshes;
    for(std::size_t submesh = 0; submesh < num_submeshes; ++submesh) {
        const std::size_t generated_partitions = partition_connected_component(graph, submesh, partition_size, min_unused_label); // this is where the magic happens
        min_unused_label += generated_partitions;
    }

    reassign_labels_by_relative_frequency(graph);
}
