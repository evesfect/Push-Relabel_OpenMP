#include "parallel.h"
#include <limits>
#include <algorithm>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

// Main parallel push-relabel algorithm
int PushRelabelParallel::maxFlow(FlowNetwork& network, int source, int sink, int num_threads) {
    int n = network.getNumVertices();
    
    std::cout << "Starting Parallel Push-Relabel with " << n << " vertices (source=" 
              << source << ", sink=" << sink << ")" << std::endl;
    
    // Ensure valid source and sink
    if (source < 0 || source >= n || sink < 0 || sink >= n || source == sink) {
        throw std::invalid_argument("Invalid source or sink");
    }
    
    // Set number of threads if specified
    #ifdef _OPENMP
    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
        std::cout << "Using " << num_threads << " threads" << std::endl;
    } else {
        std::cout << "Using " << omp_get_max_threads() << " threads" << std::endl;
    }
    #else
    std::cout << "OpenMP not enabled, running sequentially" << std::endl;
    #endif
    
    auto& graph = const_cast<std::vector<std::vector<FlowNetwork::Edge>>&>(network.getGraph());
    
    // Initialize data structures
    std::vector<int> excess(n, 0);       // Excess flow at each vertex
    std::vector<int>height(n, 0);        // Height of each vertex
    std::queue<int> active_vertices;     // Queue of active vertices
    bool* in_queue = new bool[n]();      // Track if vertex is in queue
    
    // Initialize preflow
    initialize(graph, excess, height, active_vertices, in_queue, source, sink, n);
    
    std::cout << "Initialization complete. Active vertices: " << active_vertices.size() << std::endl;
    
    // Operation counter for global relabeling
    int op_count = 0;
    int relabel_threshold = n * GLOBAL_RELABEL_FREQ;
    int iterations = 0;
    int max_iterations = MAX_ITERATIONS;
    bool done = false;
    
    // Main parallel loop
    #pragma omp parallel shared(graph, excess, height, active_vertices, in_queue, op_count, done, iterations)
    {
        // Thread-local variables
        std::vector<int> local_vertices;
        std::vector<int> local_to_add;
        
        // Reserve space to avoid reallocations
        local_vertices.reserve(BATCH_SIZE);
        local_to_add.reserve(BATCH_SIZE * 2); // May add multiple vertices per processing
        
        #pragma omp barrier // Ensure all threads are ready
        
        while (!done && iterations < max_iterations) {
            bool queue_empty = false;
            
            // Get batch of vertices from queue (critical section)
            #pragma omp critical(queue_access)
            {
                queue_empty = active_vertices.empty();
                if (!queue_empty) {
                    dequeueBatch(active_vertices, local_vertices, BATCH_SIZE);
                }
            }
            
            // Process local vertices
            if (!queue_empty && !local_vertices.empty()) {
                for (int v : local_vertices) {
                    // Skip source and sink
                    if (v != source && v != sink) {
                        discharge(graph, excess, height, local_to_add, in_queue, v, source, sink);
                    }
                }
                
                // Add any new active vertices to queue
                if (!local_to_add.empty()) {
                    #pragma omp critical(queue_access)
                    {
                        enqueueBatch(active_vertices, local_to_add, in_queue);
                    }
                    local_to_add.clear();
                }
                
                // Update operation count for global relabeling
                #pragma omp atomic
                op_count += local_vertices.size();
                
                // Clear local batch for next iteration
                local_vertices.clear();
            }
            
            // Check termination - use a reduction to check if all threads think we're done
            bool local_done = queue_empty;
            #pragma omp barrier
            
            // First, check if any thread has active vertices
            #pragma omp single
            {
                bool queue_still_empty = active_vertices.empty();
                if (!queue_still_empty) {
                    local_done = false;  // If queue has elements, we're not done
                }
            }
            
            // Make sure all threads get the same result
            #pragma omp barrier
            
            // All threads reduce to a single "done" value
            #pragma omp single
            {
                done = local_done;
                
                // Increment iteration count
                iterations++;
                
                // Print progress periodically
                if (iterations % 100 == 0) {
                    std::cout << "Iteration " << iterations << ", queue size: " 
                              << active_vertices.size() << std::endl;
                }
                
                // Check if global relabeling should be performed
                if (op_count >= relabel_threshold && !done) {
                    std::cout << "Performing global relabeling at iteration " << iterations << std::endl;
                    globalRelabeling(graph, height, excess, active_vertices, in_queue, source, sink, n);
                    op_count = 0;
                    done = active_vertices.empty();  // Recheck after relabeling
                }
            }
            
            // Make sure all threads get the updated done value
            #pragma omp barrier
        }
    }
    
    if (iterations >= max_iterations) {
        std::cout << "WARNING: Reached maximum iterations (" << max_iterations << ")" << std::endl;
    }
    
    std::cout << "Parallel Push-Relabel completed in " << iterations << " iterations" << std::endl;
    
    // Calculate the max flow
    int max_flow = 0;
    for (const auto& edge : graph[source]) {
        max_flow += edge.flow;
    }
    
    std::cout << "Max flow: " << max_flow << std::endl;
    
    // Cleanup
    delete[] in_queue;
    return max_flow;
}

// Initialize preflow by saturating edges from source
void PushRelabelParallel::initialize(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                    std::vector<int>& excess, std::vector<int>& height,
                                    std::queue<int>& active_vertices, bool* in_queue,
                                    int source, int sink, int n) {
    std::cout << "Initializing preflow..." << std::endl;
    
    // Set source height to n
    height[source] = n;
    
    // Saturate all edges from source
    for (size_t i = 0; i < graph[source].size(); i++) {
        auto& edge = graph[source][i];
        auto& rev_edge = graph[edge.to][edge.rev];
        
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

// Discharge operation - process a vertex to push its excess flow
void PushRelabelParallel::discharge(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                  std::vector<int>& excess, std::vector<int>& height,
                                  std::vector<int>& local_to_add, bool* in_queue,
                                  int u, int source, int sink) {
    // Skip source and sink
    if (u == source || u == sink) return;
    
    // Get the current excess (thread-safe read)
    int current_excess;
    #pragma omp atomic read
    current_excess = excess[u];
    
    // Safety valve - small constant to prevent infinite loops
    int max_pushes = 2 * graph[u].size();
    int pushes = 0;
    
    // Keep pushing flow until excess is gone or no admissible edge or max pushes reached
    while (current_excess > 0 && pushes < max_pushes) {
        pushes++;
        
        // Try to find an admissible edge
        bool pushed = false;
        for (size_t i = 0; i < graph[u].size(); i++) {
            auto& edge = graph[u][i];
            
            // Check if edge is admissible and has residual capacity
            int v_height;
            #pragma omp atomic read
            v_height = height[edge.to];
            
            if (edge.capacity > edge.flow && height[u] == v_height + 1) {
                // Push flow along the edge
                int push_amount = std::min(current_excess, edge.capacity - edge.flow);
                if (push_amount <= 0) continue;
                
                // Update flow values
                edge.flow += push_amount;
                graph[edge.to][edge.rev].flow -= push_amount;
                
                // Update excess
                current_excess -= push_amount;
                
                // Update the excess of the destination vertex (thread-safe)
                #pragma omp atomic
                excess[edge.to] += push_amount;
                
                // Add the destination to the active queue if needed
                if (edge.to != source && edge.to != sink && !in_queue[edge.to]) {
                    local_to_add.push_back(edge.to);
                }
                
                pushed = true;
                
                // If no more excess, break
                if (current_excess == 0) break;
            }
        }
        
        // If no push was possible, relabel the vertex
        if (!pushed) {
            // Find minimum height of neighbors with residual capacity
            int min_height = std::numeric_limits<int>::max();
            for (const auto& edge : graph[u]) {
                if (edge.capacity > edge.flow) {
                    int v_height;
                    #pragma omp atomic read
                    v_height = height[edge.to];
                    min_height = std::min(min_height, v_height);
                }
            }
            
            if (min_height != std::numeric_limits<int>::max()) {
                // Update height (thread-safe write)
                int new_height = min_height + 1;
                height[u] = new_height;
            } else {
                // No path to sink, can't push any more flow
                break;
            }
        }
    }
    
    // Update the actual excess value (thread-safe write)
    #pragma omp atomic write
    excess[u] = current_excess;
    
    // If vertex still has excess, add it back to queue
    if (current_excess > 0) {
        local_to_add.push_back(u);
    }
}

// Global relabeling using backward BFS from sink
void PushRelabelParallel::globalRelabeling(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                         std::vector<int>& height, 
                                         std::vector<int>& excess,
                                         std::queue<int>& active_vertices,
                                         bool* in_queue,
                                         int source, int sink, int n) {
    // Reset heights and visited array (sequential for simplicity)
    std::vector<int> new_heights(n, n); // Initialize all to "infinity"
    std::vector<bool> visited(n, false);
    
    // Source keeps its original height
    new_heights[source] = n;
    
    // Set sink height to 0
    new_heights[sink] = 0;
    
    // Prepare for BFS
    std::vector<int> current_level;
    std::vector<int> next_level;
    
    // Start BFS from sink
    current_level.push_back(sink);
    visited[sink] = true;
    
    int level = 0;
    
    // Process BFS levels
    while (!current_level.empty()) {
        next_level.clear();
        level++;
        
        // Process all vertices in the current level (in parallel)
        #pragma omp parallel
        {
            std::vector<int> thread_next;
            
            #pragma omp for nowait
            for (size_t i = 0; i < current_level.size(); i++) {
                int v = current_level[i];
                
                // Find all vertices that can reach v in the residual graph
                for (int u = 0; u < n; u++) {
                    if (visited[u]) continue; // Skip already visited
                    
                    for (const auto& edge : graph[u]) {
                        if (edge.to == v && edge.capacity > edge.flow) {
                            // Found path from u to v in residual graph
                            new_heights[u] = level;
                            thread_next.push_back(u);
                            visited[u] = true;
                            break;
                        }
                    }
                }
            }
            
            // Merge thread-local results safely
            #pragma omp critical(next_level)
            {
                next_level.insert(next_level.end(), thread_next.begin(), thread_next.end());
            }
        }
        
        // Move to next level
        current_level.swap(next_level);
    }
    
    // Update heights with the new values
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        height[i] = new_heights[i];
    }
    
    // Rebuild the active vertices queue
    while (!active_vertices.empty()) {
        active_vertices.pop();
    }
    
    // Reset in_queue array
    for (int i = 0; i < n; i++) {
        in_queue[i] = false;
    }
    
    // Add vertices with excess to the queue
    for (int i = 0; i < n; i++) {
        if (i != source && i != sink && excess[i] > 0) {
            active_vertices.push(i);
            in_queue[i] = true;
        }
    }
}

// Dequeue a batch of vertices from the global queue
void PushRelabelParallel::dequeueBatch(std::queue<int>& global_queue, 
                                     std::vector<int>& local_batch,
                                     int max_size) {
    int count = 0;
    while (!global_queue.empty() && count < max_size) {
        local_batch.push_back(global_queue.front());
        global_queue.pop();
        count++;
    }
}

// Add a batch of vertices to the global queue
void PushRelabelParallel::enqueueBatch(std::queue<int>& global_queue,
                                    std::vector<int>& local_batch,
                                    bool* in_queue) {
    for (int v : local_batch) {
        if (!in_queue[v]) {
            global_queue.push(v);
            in_queue[v] = true;
        }
    }
}