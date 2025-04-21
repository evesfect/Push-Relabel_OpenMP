#ifndef PARALLEL_H
#define PARALLEL_H

#include "flow_network.h"
#include <queue>
#include <vector>
#include <iostream>
#include <omp.h>
#include <boost/lockfree/queue.hpp>

class PushRelabelParallel {
public:
    // Added after experiencing infinite loops
    static const int MAX_ITERATIONS = 1000000;
    
    static int maxFlow(FlowNetwork& network, int source, int sink, int num_threads = 0);
    
private:

    static void initialize(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                          std::vector<int>& excess, std::vector<int>& height,
                          boost::lockfree::queue<int>& active_vertices, bool* in_queue,
                          int source, int sink, int n);
    
    static bool push(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                    std::vector<int>& excess, boost::lockfree::queue<int>& active_vertices,
                    bool* in_queue, int u, int v_idx);
    
    static bool relabel(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                       std::vector<int>& height, int u);
    
    static void discharge(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                        std::vector<int>& excess, std::vector<int>& height,
                        boost::lockfree::queue<int>& active_vertices, bool* in_queue,
                        int u, int source, int sink, int& total_relabels);
    
};

#endif