#ifndef SEQUENTIAL_H
#define SEQUENTIAL_H

#include "flow_network.h"
#include <vector>
#include <iostream>
#include <queue>
#include <limits>

class PushRelabelSequential {
public:
    static const int MAX_ITERATIONS = 100000000;
    static const int INF_HEIGHT = std::numeric_limits<int>::max(); // For GL

    static int maxFlow(FlowNetwork& network, int source, int sink);
    static void globalRelabel(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                             std::vector<int>& height, int source, int sink, int n);

private:
    static void initialize(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                          std::vector<int>& excess, std::vector<int>& height,
                          std::queue<int>& active_vertices, std::vector<bool>& in_queue,
                          int source, int sink, int n);

    static bool push(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                    std::vector<int>& excess, std::vector<int>& height,
                    std::queue<int>& active_vertices, std::vector<bool>& in_queue,
                    int u, int v_idx, int source, int sink);

    static bool relabel(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                       std::vector<int>& height, int u, int n);

    static int discharge(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                        std::vector<int>& excess, std::vector<int>& height,
                        std::queue<int>& active_vertices, std::vector<bool>& in_queue,
                        int u, int source, int sink, int n);

    static void printState(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                          const std::vector<int>& excess, const std::vector<int>& height,
                          const std::queue<int>& active_vertices, int n);
};

#endif