#ifndef GENERATOR_H
#define GENERATOR_H

#include "flow_network.h"
#include <random>
#include <string>

class FlowNetworkGenerator {
public:
    FlowNetworkGenerator(unsigned int seed = std::random_device{}());

    // Generate random flow network
    FlowNetwork generateRandom(int num_vertices, double edge_probability, 
                             int min_capacity, int max_capacity, 
                             int source = 0, int sink = -1);

    // Generate layered network (useful for testing)
    FlowNetwork generateLayered(int num_layers, int vertices_per_layer,
                              double edge_probability,
                              int min_capacity, int max_capacity);
                              
    // Generate a dense network
    FlowNetwork generateDense(int num_vertices, int min_capacity, int max_capacity,
                            int source = 0, int sink = -1);

private:
    std::mt19937 rng;
};

#endif // GENERATOR_H