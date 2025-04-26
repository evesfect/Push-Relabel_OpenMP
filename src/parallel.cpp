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
    
    std::vector<int> excess(n, 0);               // Current excess flow of each vertex
    std::vector<int> height(n, 0);               // Current height of each vertex
    std::vector<int> excess_changes(n, 0);       // Temporary storage for excess changes
    std::vector<int> new_heights(n, 0);          // Temporary storage for new heights
    
    std::vector<int> active_vertices(n);         // Current active vertices
    std::vector<bool> is_active(n, false);       // Track if vertex is active
    std::vector<bool> still_active(n, false);    // Track if vertex is still active after push
    std::vector<int> next_active_vertices(n);    // Active vertices for next iteration
    int active_count = 0;                        // Number of active vertices
    int next_active_count = 0;                   // Number of active vertices for next iteration
    
    initialize(graph, excess, height, active_vertices, is_active, active_count, source, sink, n);
    
    int iterations = 0;
    
    // Main loop - continue until no active vertices remain
    while (active_count > 0 && iterations < MAX_ITERATIONS) {
        // Clear temporary structures
        std::fill(excess_changes.begin(), excess_changes.end(), 0);
        std::fill(new_heights.begin(), new_heights.end(), 0);
        std::fill(still_active.begin(), still_active.end(), false);
        next_active_count = 0;
        
        // Phase 1: Process all active vertices and compute excess changes
        pushPhase(graph, excess, height, excess_changes, active_vertices, still_active, active_count, source, sink);
        
        // Phase 2: Compute new labels for vertices that are still active
        labelComputationPhase(graph, excess, height, new_heights, still_active);
        
        // Phase 3: Apply the new labels
        labelApplicationPhase(height, new_heights, still_active);
        
        // Phase 4: Apply excess changes and update active vertices for next iteration
        excessUpdatePhase(excess, excess_changes, next_active_vertices, is_active, next_active_count, source, sink, n);
        
        // Swap active vertices for next iteration
        active_vertices.swap(next_active_vertices);
        active_count = next_active_count;
        
        iterations++;
    }
    
    if (iterations >= MAX_ITERATIONS) {
        std::cout << "WARNING: Reached maximum iterations (" << MAX_ITERATIONS << ")." << std::endl;
    }
    
    // Calculate the max flow 
    int max_flow = 0;
    for (const auto& edge : graph[source]) {
        max_flow += edge.flow;
    }
    
    std::cout << "Completed in " << iterations << " iterations" << std::endl;
    std::cout << "Max flow: " << max_flow << std::endl;
    
    return max_flow;
}

// Initialize preflow - saturate edges from source
void PushRelabelParallel::initialize(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                     std::vector<int>& excess, std::vector<int>& height,
                                     std::vector<int>& active_vertices, std::vector<bool>& is_active,
                                     int& active_count, int source, int sink, int n) {
    // Set source height to n
    height[source] = n;
    active_count = 0;
    
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
        
        // Add to active vertices if it's not source or sink and has excess
        if (edge.to != source && edge.to != sink && excess[edge.to] > 0 && !is_active[edge.to]) {
            is_active[edge.to] = true;
            active_vertices[active_count++] = edge.to;
        }
    }
}

// Phase 1: Push operations - store changes in temporary arrays
void PushRelabelParallel::pushPhase(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                    const std::vector<int>& excess, const std::vector<int>& height,
                                    std::vector<int>& excess_changes, const std::vector<int>& active_vertices,
                                    std::vector<bool>& still_active, int active_count, int source, int sink) {
    // Process each active vertex
    for (int i = 0; i < active_count; i++) {
        int u = active_vertices[i];
        if (u == source || u == sink) continue;
        
        bool pushed = false;
        
        // Check all edges for possible pushes
        for (size_t j = 0; j < graph[u].size(); j++) {
            auto& edge = graph[u][j];
            
            // Check if edge is admissible
            if (edge.capacity > edge.flow && height[u] == height[edge.to] + 1) {
                // Attempt to push flow
                if (push(graph, excess, excess_changes, u, j)) {
                    pushed = true;
                }
            }
        }
        
        // If vertex is still active after all push attempts
        if (excess[u] + excess_changes[u] > 0) {
            still_active[u] = true;
        }
    }
}

// Push operation - updates edge flows and records excess changes
bool PushRelabelParallel::push(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                              const std::vector<int>& excess, std::vector<int>& excess_changes,
                              int u, int v_idx) {
    auto& edge = graph[u][v_idx];
    auto& rev_edge = graph[edge.to][edge.rev];
    
    // Calculate flow to push (minimum of excess and residual capacity)
    int flow = std::min(excess[u] + excess_changes[u], edge.capacity - edge.flow);
    
    if (flow <= 0) {
        return false;  // No flow was pushed
    }
    
    // Update flow values
    edge.flow += flow;
    rev_edge.flow -= flow;
    
    // Update excess changes
    excess_changes[u] -= flow;
    excess_changes[edge.to] += flow;
    
    return true;  // Flow was pushed successfully
}

// Phase 2: Compute new labels for vertices that are still active
void PushRelabelParallel::labelComputationPhase(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                              const std::vector<int>& excess, const std::vector<int>& height,
                                              std::vector<int>& new_heights, const std::vector<bool>& still_active) {
    // Compute new heights for vertices that are still active
    for (int u = 0; u < still_active.size(); u++) {
        if (still_active[u]) {
            int new_height;
            if (computeNewLabel(graph, height, u, new_height)) {
                new_heights[u] = new_height;
            } else {
                // If relabeling fails, keep the old height
                new_heights[u] = height[u];
            }
        }
    }
}

// Compute new label for a vertex
bool PushRelabelParallel::computeNewLabel(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                        const std::vector<int>& height, int u, int& new_height) {
    int min_height = std::numeric_limits<int>::max();
    
    // Find minimum height of neighboring vertices with residual capacity
    for (const auto& edge : graph[u]) {
        if (edge.capacity > edge.flow) {
            min_height = std::min(min_height, height[edge.to]);
        }
    }
    
    // Update height to be one more than minimum
    if (min_height != std::numeric_limits<int>::max()) {
        new_height = min_height + 1;
        return true;
    }
    
    return false;  // Could not relabel
}

// Phase 3: Apply the new labels
void PushRelabelParallel::labelApplicationPhase(std::vector<int>& height, const std::vector<int>& new_heights,
                                              const std::vector<bool>& still_active) {
    // Apply new heights to vertices that are still active
    for (int u = 0; u < still_active.size(); u++) {
        if (still_active[u]) {
            if (new_heights[u] > height[u]) {
                height[u] = new_heights[u];
            }
        }
    }
}

// Phase 4: Apply excess changes and update active vertices
void PushRelabelParallel::excessUpdatePhase(std::vector<int>& excess, const std::vector<int>& excess_changes,
                                          std::vector<int>& next_active_vertices, std::vector<bool>& is_active,
                                          int& next_active_count, int source, int sink, int n) {
    // Reset active tracking
    std::fill(is_active.begin(), is_active.end(), false);
    next_active_count = 0;
    
    // Apply excess changes and determine active vertices for next iteration
    for (int i = 0; i < n; i++) {
        // Skip source and sink
        if (i == source || i == sink) continue;
        
        // Apply excess change
        excess[i] += excess_changes[i];
        
        // If vertex has excess, add to next active set
        if (excess[i] > 0) {
            is_active[i] = true;
            next_active_vertices[next_active_count++] = i;
        }
    }
}

// Debug method to print current state
void PushRelabelParallel::printState(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                    const std::vector<int>& excess, const std::vector<int>& height,
                                    const std::vector<int>& active_vertices, int active_count, int n) {
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
    
    std::cout << "  Active vertices: ";
    for (int i = 0; i < active_count; i++) {
        std::cout << active_vertices[i] << " ";
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