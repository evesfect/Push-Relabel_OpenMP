#include "generator.h"

FlowNetworkGenerator::FlowNetworkGenerator(unsigned int seed) : rng(seed) {}

FlowNetwork FlowNetworkGenerator::generateRandom(int num_vertices, double edge_probability, 
                                         int min_capacity, int max_capacity,
                                         int source, int sink) {
    if (sink == -1) sink = num_vertices - 1;
    
    FlowNetwork network(num_vertices);
    std::uniform_real_distribution<double> edge_dist(0.0, 1.0);
    std::uniform_int_distribution<int> capacity_dist(min_capacity, max_capacity);

    for (int from = 0; from < num_vertices; ++from) {
        for (int to = 0; to < num_vertices; ++to) {
            if (from != to && from < to && edge_dist(rng) < edge_probability) {
                network.addEdge(from, to, capacity_dist(rng));
            }
        }
    }

    return network;
}
