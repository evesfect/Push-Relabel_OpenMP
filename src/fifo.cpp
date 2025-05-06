#include "fifo.h"
#include <limits>
#include <algorithm>
#include <iomanip>


int FIFO::maxFlow(FlowNetwork& network, int source, int sink) {
    int n = network.getNumVertices();
    
    std::cout << "Starting with " << n << " vertices (source=" 
              << source << ", sink=" << sink << ")" << std::endl;
    
    if (source < 0 || source >= n || sink < 0 || sink >= n || source == sink) {
        throw std::invalid_argument("Invalid source or sink");
    }
    
    auto& graph = const_cast<std::vector<std::vector<FlowNetwork::Edge>>&>(network.getGraph());
    
    std::vector<int> excess(n, 0);
    std::vector<int> height(n, 0);
    std::queue<int> active_vertices;
    bool* in_queue = new bool[n]();
    
    initialize(graph, excess, height, active_vertices, in_queue, source, sink, n);
    
    int iterations = 0;
    int total_relabels = 0;
    
    while (!active_vertices.empty() && iterations < MAX_ITERATIONS) {
        //simple implementation*
        int u = active_vertices.front();
        
        active_vertices.pop();
        in_queue[u] = false;
        
        discharge(graph, excess, height, active_vertices, in_queue, u, source, sink, total_relabels);
        
        iterations++;
    }
    
    if (iterations >= MAX_ITERATIONS) {
        std::cout << "WARNING: Reached maximum iterations (" << MAX_ITERATIONS 
                  << ")." << std::endl;
    }
    
    int max_flow = 0;
    for (const auto& edge : graph[source]) {
        max_flow += edge.flow;
    }
    
    std::cout << "Completed in " << iterations << " iterations" << std::endl;
    std::cout << "Max flow: " << max_flow << std::endl;
    
    // Sanity check
    /*
    int sink_flow = 0;
    for (size_t i = 0; i < graph.size(); i++) {
        for (const auto& edge : graph[i]) {
            if (edge.to == sink && edge.flow > 0) {
                sink_flow += edge.flow;
            }
        }
    }
    std::cout << "Flow into sink: " << sink_flow << " (should match max flow)" << std::endl;
    */
    delete[] in_queue;
    return max_flow;
}


void FIFO::initialize(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                     std::vector<int>& excess, std::vector<int>& height,
                                     std::queue<int>& active_vertices, bool* in_queue,
                                     int source, int sink, int n) {
    
    height[source] = n;
    
    for (size_t i = 0; i < graph[source].size(); i++) {
        auto& edge = const_cast<FlowNetwork::Edge&>(graph[source][i]);
        auto& rev_edge = const_cast<FlowNetwork::Edge&>(graph[edge.to][edge.rev]);
        
        int flow = edge.capacity;
        edge.flow = flow;
        rev_edge.flow = -flow;
        
        excess[edge.to] += flow;

        if (edge.to != source && edge.to != sink && excess[edge.to] > 0 && !in_queue[edge.to]) {
            active_vertices.push(edge.to);
            in_queue[edge.to] = true;
        }
    }
    
    //std::cout << "Init done. Active vertices: " << active_vertices.size() << std::endl;
}

bool FIFO::push(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                               std::vector<int>& excess, std::queue<int>& active_vertices,
                               bool* in_queue, int u, int v_idx) {
    auto& edge = graph[u][v_idx];
    auto& rev_edge = graph[edge.to][edge.rev];
    
    int flow = std::min(excess[u], edge.capacity - edge.flow);
    
    if (flow <= 0) {
        return false;
    }
    
    edge.flow += flow;
    rev_edge.flow -= flow;
    
    excess[u] -= flow;
    excess[edge.to] += flow;
    
    if (excess[edge.to] > 0 && !in_queue[edge.to]) {
        active_vertices.push(edge.to);
        in_queue[edge.to] = true;
    }
    
    return true;
}

bool FIFO::relabel(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                  std::vector<int>& height, int u) {
    int old_height = height[u];
    int min_height = std::numeric_limits<int>::max();
    
    for (const auto& edge : graph[u]) {
        if (edge.capacity > edge.flow) {
            min_height = std::min(min_height, height[edge.to]);
        }
    }
    
    if (min_height != std::numeric_limits<int>::max()) {
        height[u] = min_height + 1;
        return true;
    }
    
    std::cout << "WARNING: Could not relabel vertex " << u 
              << " - no neighbors with residual capacity" << std::endl;
    
    return false;
}

void FIFO::discharge(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                    std::vector<int>& excess, std::vector<int>& height,
                                    std::queue<int>& active_vertices, bool* in_queue,
                                    int u, int source, int sink, int& total_relabels) {
    if (u == source || u == sink) return;
    
    int discharge_iterations = 0;
    const int MAX_DISCHARGE_ITERATIONS = 2 * graph.size();
    
    while (excess[u] > 0 && discharge_iterations < MAX_DISCHARGE_ITERATIONS) {
        bool pushed = false;
        for (size_t i = 0; i < graph[u].size(); i++) {
            auto& edge = graph[u][i];
            
            if (edge.capacity > edge.flow && height[u] == height[edge.to] + 1) {
                pushed = push(graph, excess, active_vertices, in_queue, u, i);
                if (pushed) break;
            }
        }
        
        // If cannot find edge, relabel
        // Error check is obsolete here, but had issues with relabeling previously
        // Kept for debugging
        if (!pushed) {
            bool relabeled = relabel(graph, height, u);
            if (relabeled) {
                total_relabels++;
            } else {

                std::cout << "ERROR: Could not relabel vertex " << u 
                          << " with excess " << excess[u] << std::endl;
                break;
            }
        }
        
        discharge_iterations++;
    }
    
    if (discharge_iterations >= MAX_DISCHARGE_ITERATIONS) {
        std::cout << "WARNING: Reached maximum discharge iterations for vertex " 
                  << u << ", excess=" << excess[u] << std::endl;
    }
    
    // If vertex has excess after discharge, add back to the queue
    if (excess[u] > 0 && u != source && u != sink && !in_queue[u]) {
        active_vertices.push(u);
        in_queue[u] = true;
    }
}

void FIFO::printState(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                     const std::vector<int>& excess, const std::vector<int>& height,
                                     const std::queue<int>& active_vertices, int n) {
    std::cout << "Current State:" << std::endl;

    std::cout << "  Heights: ";
    for (int i = 0; i < n; i++) {
        std::cout << height[i] << " ";
    }
    std::cout << std::endl;
    
    std::cout << "  Excess: ";
    for (int i = 0; i < n; i++) {
        std::cout << excess[i] << " ";
    }
    std::cout << std::endl;
    
    std::cout << "  Flow network:" << std::endl;
    for (int u = 0; u < n; u++) {
        for (const auto& edge : graph[u]) {
            if (edge.capacity > 0) {  // Only print forward edges
                std::cout << "    " << u << " -> " << edge.to 
                          << " (flow=" << edge.flow << "/" << edge.capacity << ")" << std::endl;
            }
        }
    }
    std::cout << std::endl;
}