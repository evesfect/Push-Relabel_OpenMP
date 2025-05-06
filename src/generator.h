#ifndef GENERATOR_H
#define GENERATOR_H

#include "flow_network.h"
#include <random>
#include <string>

class FlowNetworkGenerator {
public:
    FlowNetworkGenerator(unsigned int seed = std::random_device{}());

    FlowNetwork generateRandom(int num_vertices, double edge_probability, 
                             int min_capacity, int max_capacity, 
                             int source = 0, int sink = -1);


private:
    std::mt19937 rng;
};

#endif