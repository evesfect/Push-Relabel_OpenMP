#include "parallel.h"
#include <limits>
#include <algorithm>
#include <iomanip>

#ifdef _OPENMP
#include <omp.h>
#endif

int PushRelabelParallel::maxFlow(FlowNetwork& network, int source, int sink, int num_threads) {
    int n = network.getNumVertices();
    
    std::cout << "Starting parallel push-relabel with " << n << " vertices (source=" 
              << source << ", sink=" << sink << ")" << std::endl;
    
    // Set number of threads if specified
    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
    }
    
    // Validate source and sink
    if (source < 0 || source >= n || sink < 0 || sink >= n || source == sink) {
        throw std::invalid_argument("Invalid source or sink");
    }
    
    auto& graph = const_cast<std::vector<std::vector<FlowNetwork::Edge>>&>(network.getGraph());
    
    std::vector<int> excess(n, 0);
    std::vector<int> height(n, 0);
    boost::lockfree::queue<int> active_vertices(n);
    bool* in_queue = new bool[n]();
    
    // Initialize preflow
    initialize(graph, excess, height, active_vertices, in_queue, source, sink, n);
    
    int iterations = 0;
    int total_relabels = 0;
    
    // Main Loop
    // Works on vertices with excess flow
    int u;
    while (active_vertices.pop(u) && iterations < MAX_ITERATIONS) {
        // Choose an active vertex
        in_queue[u] = false;
        
        discharge(graph, excess, height, active_vertices, in_queue, u, source, sink, total_relabels);
        
        iterations++;
    }
    
    if (iterations >= MAX_ITERATIONS) {
        std::cout << "WARNING: Reached maximum iterations (" << MAX_ITERATIONS 
                  << ")." << std::endl;
    }
    
    // Calculate the max flow 
    int max_flow = 0;
    for (const auto& edge : graph[source]) {
        max_flow += edge.flow;
    }
    
    std::cout << "Completed in " << iterations << " iterations" << std::endl;
    std::cout << "Max flow: " << max_flow << std::endl;
    
    delete[] in_queue;
    return max_flow;
}

// Initialize preflow - saturate edges from source
void PushRelabelParallel::initialize(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                     std::vector<int>& excess, std::vector<int>& height,
                                     boost::lockfree::queue<int>& active_vertices, bool* in_queue,
                                     int source, int sink, int n) {
    
    // Set source height to n
    height[source] = n;
    
    // Saturate all edges from source
    for (size_t i = 0; i < graph[source].size(); i++) {
        auto& edge = const_cast<FlowNetwork::Edge&>(graph[source][i]);
        auto& rev_edge = const_cast<FlowNetwork::Edge&>(graph[edge.to][edge.rev]);
        
        // Push maximum flow from source to neighbor
        int flow = edge.capacity;
        edge.flow = flow;
        rev_edge.flow = -flow;  // Update reverse edge
        
        // Update excess flow
        excess[edge.to] += flow;
        
        // Add neighbor to queue if it's not source or sink and has excess
        if (edge.to != source && edge.to != sink && excess[edge.to] > 0 && !in_queue[edge.to]) {
            active_vertices.push(edge.to);
            in_queue[edge.to] = true;
        }
    }
}

// Push operation - push flow from u to v
bool PushRelabelParallel::push(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                               std::vector<int>& excess, boost::lockfree::queue<int>& active_vertices,
                               bool* in_queue, int u, int v_idx) {
    auto& edge = graph[u][v_idx];
    auto& rev_edge = graph[edge.to][edge.rev];
    
    // Calculate flow to push (minimum of excess and residual capacity)
    int flow = std::min(excess[u], edge.capacity - edge.flow);
    
    if (flow <= 0) {
        return false;  // No flow was pushed
    }
    
    // Update flow values
    edge.flow += flow;
    rev_edge.flow -= flow;
    
    // Update excess flow
    excess[u] -= flow;
    excess[edge.to] += flow;
    
    // Add vertex to queue if it has excess and is not already in queue
    if (excess[edge.to] > 0 && !in_queue[edge.to]) {
        active_vertices.push(edge.to);
        in_queue[edge.to] = true;
    }
    
    return true;  // Flow was pushed successfully
}

// Relabel operation - increase height of vertex u
bool PushRelabelParallel::relabel(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                  std::vector<int>& height, int u) {
    int old_height = height[u];
    int min_height = std::numeric_limits<int>::max();
    
    // Find minimum height of neighboring vertices with residual capacity
    for (const auto& edge : graph[u]) {
        if (edge.capacity > edge.flow) {
            min_height = std::min(min_height, height[edge.to]);
        }
    }
    
    // Update height to be one more than minimum
    if (min_height != std::numeric_limits<int>::max()) {
        height[u] = min_height + 1;
        return true;  // Relabel was successful
    }
    
    std::cout << "WARNING: Could not relabel vertex " << u 
              << " - no neighbors with residual capacity" << std::endl;
    
    return false;  // Could not relabel
}

void PushRelabelParallel::discharge(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                    std::vector<int>& excess, std::vector<int>& height,
                                    boost::lockfree::queue<int>& active_vertices, bool* in_queue,
                                    int u, int source, int sink, int& total_relabels) {
    // Skip the source and the sink
    if (u == source || u == sink) return;
    
    // For infinite loops
    int discharge_iterations = 0;
    const int MAX_DISCHARGE_ITERATIONS = 2 * graph.size();
    
    // Loop until vertex has no more excess
    while (excess[u] > 0 && discharge_iterations < MAX_DISCHARGE_ITERATIONS) {
        // Find an edge suitable for push (compare heights)
        bool pushed = false;
        for (size_t i = 0; i < graph[u].size(); i++) {
            auto& edge = graph[u][i];
            
            if (edge.capacity > edge.flow && height[u] == height[edge.to] + 1) {
                pushed = push(graph, excess, active_vertices, in_queue, u, i);
                if (pushed) break; // Found edge, pushed flow
            }
        }
        
        // If cannot find edge, relabel
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