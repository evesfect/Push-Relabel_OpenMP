#include "flow_network.h"
#include <queue>
#include <vector>
#include <iostream>

class FIFO {
public:
    // Added after experiencing infinite loops
    // FIFO exceeds this limit easily around 10k nodes
    static const int MAX_ITERATIONS = 100000000;
    
    static int maxFlow(FlowNetwork& network, int source, int sink);
    
private:
    static void initialize(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                          std::vector<int>& excess, std::vector<int>& height,
                          std::queue<int>& active_vertices, bool* in_queue,
                          int source, int sink, int n);
    
    static bool push(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                    std::vector<int>& excess, std::queue<int>& active_vertices,
                    bool* in_queue, int u, int v_idx);
    
    static bool relabel(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                       std::vector<int>& height, int u);
    
    static void discharge(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                        std::vector<int>& excess, std::vector<int>& height,
                        std::queue<int>& active_vertices, bool* in_queue,
                        int u, int source, int sink, int& total_relabels);
    
    static void printState(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                          const std::vector<int>& excess, const std::vector<int>& height,
                          const std::queue<int>& active_vertices, int n);
};
