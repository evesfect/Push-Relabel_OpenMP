
#include "flow_network.h"
#include <queue>
#include <vector>
#include <iostream>

class FIFO {
public:
    // Added after experiencing infinite loops
    static const int MAX_ITERATIONS = 100000000;
    
    /// @brief Calculates maximum possible flow of a given flow network
    /// @param network 
    /// @param source 
    /// @param sink 
    /// @return 
    static int maxFlow(FlowNetwork& network, int source, int sink);
    
private:
    /// @brief Initializes the preflow 
    /// @param graph 
    /// @param excess 
    /// @param height 
    /// @param active_vertices
    /// @param in_queue 
    /// @param source 
    /// @param sink 
    /// @param n
    static void initialize(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                          std::vector<int>& excess, std::vector<int>& height,
                          std::queue<int>& active_vertices, bool* in_queue,
                          int source, int sink, int n);
    
    /// @brief Pushes flow from u to v
    /// @param graph 
    /// @param excess 
    /// @param active_vertices 
    /// @param in_queue 
    /// @param u 
    /// @param v_idx 
    /// @return 
    static bool push(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                    std::vector<int>& excess, std::queue<int>& active_vertices,
                    bool* in_queue, int u, int v_idx);
    
    /// @brief Relabels (increases height) of vertex u
    /// @param graph 
    /// @param height 
    /// @param u 
    /// @return 
    static bool relabel(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                       std::vector<int>& height, int u);
    
    /// @brief Pushes (discharges) excess flow from vertex u
    /// @param graph 
    /// @param excess 
    /// @param height 
    /// @param active_vertices 
    /// @param in_queue 
    /// @param u 
    /// @param source 
    /// @param sink 
    /// @param total_relabels 
    static void discharge(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                        std::vector<int>& excess, std::vector<int>& height,
                        std::queue<int>& active_vertices, bool* in_queue,
                        int u, int source, int sink, int& total_relabels);
    
    /// @brief Basic debugging tool
    /// @param graph 
    /// @param excess 
    /// @param height 
    /// @param active_vertices 
    /// @param n 
    static void printState(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                          const std::vector<int>& excess, const std::vector<int>& height,
                          const std::queue<int>& active_vertices, int n);
};
