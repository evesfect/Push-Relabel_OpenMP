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
            // Avoid self-loops and ensure forward flow (directed acyclic for simplicity)
            if (from != to && from < to && edge_dist(rng) < edge_probability) {
                network.addEdge(from, to, capacity_dist(rng));
            }
        }
    }

    return network;
}

FlowNetwork FlowNetworkGenerator::generateLayered(int num_layers, int vertices_per_layer,
                                         double edge_probability,
                                         int min_capacity, int max_capacity) {
    int num_vertices = num_layers * vertices_per_layer;
    FlowNetwork network(num_vertices);
    
    std::uniform_real_distribution<double> edge_dist(0.0, 1.0);
    std::uniform_int_distribution<int> capacity_dist(min_capacity, max_capacity);

    // Connect vertices between consecutive layers
    for (int layer = 0; layer < num_layers - 1; ++layer) {
        int start_current = layer * vertices_per_layer;
        int start_next = (layer + 1) * vertices_per_layer;
        
        for (int i = 0; i < vertices_per_layer; ++i) {
            int from = start_current + i;
            
            for (int j = 0; j < vertices_per_layer; ++j) {
                int to = start_next + j;
                
                if (edge_dist(rng) < edge_probability) {
                    network.addEdge(from, to, capacity_dist(rng));
                }
            }
        }
    }

    return network;
}

FlowNetwork FlowNetworkGenerator::generateDense(int num_vertices, int min_capacity, int max_capacity,
                                       int source, int sink) {
    if (sink == -1) sink = num_vertices - 1;
    
    FlowNetwork network(num_vertices);
    std::uniform_int_distribution<int> capacity_dist(min_capacity, max_capacity);

    // Create a dense network where almost all possible edges exist
    for (int from = 0; from < num_vertices; ++from) {
        for (int to = 0; to < num_vertices; ++to) {
            // Avoid self-loops
            if (from != to) {
                network.addEdge(from, to, capacity_dist(rng));
            }
        }
    }

    return network;
}