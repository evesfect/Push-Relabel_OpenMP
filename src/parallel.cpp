#include "parallel.h"
#include <limits>
#include <algorithm>
#include <iomanip>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

int PushRelabelParallel::maxFlow(FlowNetwork& network, int source, int sink, int num_threads) {
    int n = network.getNumVertices();
    std::cout << "Starting Parallel Push-Relabel with " << n
              << " vertices (source=" << source << ", sink=" << sink << ")" << std::endl;

    // Thread setup
    if (num_threads > 0) {
        #ifdef _OPENMP
        omp_set_num_threads(num_threads);
        #else
        std::cout << "Warning: OpenMP not enabled, running sequentially." << std::endl;
        #endif
    } else {
        #ifdef _OPENMP
        num_threads = omp_get_max_threads();
        omp_set_num_threads(num_threads);
        std::cout << "Using default OpenMP threads: " << num_threads << std::endl;
        #else
        num_threads = 1;
        std::cout << "Warning: OpenMP not enabled, running sequentially." << std::endl;
        #endif
    }

    if (source < 0 || source >= n || sink < 0 || sink >= n || source == sink) {
        throw std::invalid_argument("Invalid source or sink");
    }

    auto& graph = const_cast<std::vector<std::vector<FlowNetwork::Edge>>&>(network.getGraph());
    std::vector<int> excess(n, 0), height(n, 0);
    std::queue<int> active_vertices;
    std::vector<bool> in_queue(n, false);

    // Global relabel frequency setup
    int relabel_since_last = 0;
    const int global_freq = n;
    int total_global = 0;

    initialize(graph, excess, height, active_vertices, in_queue, source, sink, n);

    // Initial global relabel
    globalRelabel(graph, height, source, sink, n);
    total_global++;
    active_vertices = {};
    std::fill(in_queue.begin(), in_queue.end(), false);
    for (int i = 0; i < n; ++i) {
        if (i != source && i != sink && excess[i] > 0 && height[i] < n) {
            active_vertices.push(i);
            in_queue[i] = true;
        }
    }

    int iterations = 0, total_relabels = 0;
    while (!active_vertices.empty() && iterations < MAX_ITERATIONS) {
        if (relabel_since_last >= global_freq) {
            // Global relabel
            globalRelabel(graph, height, source, sink, n);
            total_global++;
            relabel_since_last = 0;

            std::queue<int> new_queue;
            std::fill(in_queue.begin(), in_queue.end(), false);
            while (!active_vertices.empty()) {
                int u = active_vertices.front(); active_vertices.pop();
                if (u != source && u != sink && excess[u] > 0 && height[u] < n) {
                    new_queue.push(u);
                    in_queue[u] = true;
                }
            }
            for (int i = 0; i < n; ++i) {
                if (i != source && i != sink && excess[i] > 0 && height[i] < n && !in_queue[i]) {
                    new_queue.push(i);
                    in_queue[i] = true;
                }
            }
            active_vertices = std::move(new_queue);
            if (active_vertices.empty()) break;
        }

        int u = active_vertices.front();
        active_vertices.pop();
        in_queue[u] = false;

        int relabels = discharge(graph, excess, height, active_vertices, in_queue, u, source, sink, n);
        total_relabels += relabels;
        relabel_since_last += relabels;
        iterations++;
    }

    if (iterations >= MAX_ITERATIONS) {
        std::cout << "WARNING: Reached maximum iterations. Max flow might be incorrect." << std::endl;
    }

    int max_flow = excess[sink];
    std::cout << "Completed in " << iterations << " iterations. Total local relabels: "
              << total_relabels << ", global relabels: " << total_global
              << ". Max flow: " << max_flow << std::endl;
    return max_flow;
}

// Initialize preflow
void PushRelabelParallel::initialize(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                     std::vector<int>& excess, std::vector<int>& height,
                                     std::queue<int>& active_vertices, std::vector<bool>& in_queue,
                                     int source, int sink, int n) {
    height.assign(n, 0);
    excess.assign(n, 0);
    in_queue.assign(n, false);
    height[source] = n;

    for (auto& edge : const_cast<std::vector<FlowNetwork::Edge>&>(graph[source])) {
        if (edge.capacity > 0) {
            int flow = edge.capacity;
            int v = edge.to;
            auto& rev = const_cast<std::vector<FlowNetwork::Edge>&>(graph[v])[edge.rev];
            edge.flow += flow;
            rev.flow -= flow;
            excess[v] += flow;
            excess[source] -= flow;
            if (v != source && v != sink && !in_queue[v]) {
                active_vertices.push(v);
                in_queue[v] = true;
            }
        }
    }
}

// Push operation (single edge)
bool PushRelabelParallel::push(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                               std::vector<int>& excess, std::vector<int>& height,
                               std::queue<int>& active_vertices, std::vector<bool>& in_queue,
                               int u, int idx, int source, int sink) {
    auto& e = graph[u][idx];
    int v = e.to;
    if (e.capacity - e.flow <= 0 || height[u] != height[v] + 1) return false;
    auto& rev = graph[v][e.rev];
    int delta = std::min(excess[u], e.capacity - e.flow);
    if (delta <= 0) return false;
    e.flow += delta;
    rev.flow -= delta;
    excess[u] -= delta;
    excess[v] += delta;
    if (v != source && v != sink && excess[v] > 0 && !in_queue[v]) {
        active_vertices.push(v);
        in_queue[v] = true;
    }
    return true;
}

// Relabel node u
bool PushRelabelParallel::relabel(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                  std::vector<int>& height, int u, int n) {
    int min_h = INF_HEIGHT;
    for (auto& e : graph[u]) {
        if (e.capacity - e.flow > 0)
            min_h = std::min(min_h, height[e.to]);
    }
    height[u] = (min_h < INF_HEIGHT ? min_h + 1 : n);
    return true;
}

// Discharge excess at u
int PushRelabelParallel::discharge(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                   std::vector<int>& excess, std::vector<int>& height,
                                   std::queue<int>& active_vertices, std::vector<bool>& in_queue,
                                   int u, int source, int sink, int n) {
    if (u == source || u == sink || excess[u] <= 0) return 0;
    int relabels = 0;
    while (excess[u] > 0 && height[u] < n) {
        bool did_push = false;
        for (size_t i = 0; i < graph[u].size() && excess[u] > 0; ++i) {
            if (push(graph, excess, height, active_vertices, in_queue, u, i, source, sink))
                did_push = true;
        }
        if (!did_push && excess[u] > 0) {
            relabel(graph, height, u, n);
            relabels++;
        }
    }
    return relabels;
}

// Global relabel via backward BFS
void PushRelabelParallel::globalRelabel(std::vector<std::vector<FlowNetwork::Edge>>& graph, // Needs non-const graph
                                         std::vector<int>& height, int source, int sink, int n) {

    // --- Start Timing (Optional) ---
    // double start_time = 0;
    // #ifdef _OPENMP
    // start_time = omp_get_wtime();
    // #endif
    // --- End Timing ---

    // 1. Initialization (Parallel)
    // Use n as the initial "infinity" / unvisited marker.
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        height[i] = n;
    }
    // Set sink height and ensure it's marked correctly if n=0
    if (n > 0) {
       height[sink] = 0;
    }


    // 2. Level-Synchronous BFS Setup
    std::vector<int> frontier;
    if (n > 0) { // Avoid pushing sink if n=0
       frontier.push_back(sink);
    }
    frontier.reserve(n); // Reserve space, helps avoid reallocations

    std::vector<int> next_frontier;
    next_frontier.reserve(n); // Reserve space

    int current_level = 0;
    int num_threads = 1; // Default to 1 if OpenMP is not used
    #ifdef _OPENMP
    num_threads = omp_get_max_threads(); // Get number of threads available
    #endif
    std::vector<std::vector<int>> local_next_frontiers(num_threads); // Thread-local lists

    // 3. BFS Loop: Process level by level
    while (!frontier.empty()) {
        // Clear thread-local lists for the new level
        for (int i = 0; i < num_threads; ++i) {
            local_next_frontiers[i].clear();
            // Optional: Reserve estimated space in local lists if possible
        }
        next_frontier.clear(); // Clear the global next frontier from previous iteration

        #pragma omp parallel // Start parallel region for processing the current frontier
        {
            int tid = 0;
            #ifdef _OPENMP
            tid = omp_get_thread_num();
            #endif

            // Parallel loop over the current frontier nodes (v)
            // Dynamic scheduling can be beneficial if processing time per node varies
            #pragma omp for schedule(dynamic, 128) // Chunk size 128 as a starting point
            for (size_t i = 0; i < frontier.size(); ++i) {
                 int v = frontier[i];

                 // Explore edges *entering* v in the residual graph.
                 // Iterate through neighbors 'u' of 'v' (edge v->u exists).
                 // Check the residual capacity of the *reverse* edge (u->v).
                 for (const auto& edge_from_v : graph[v]) { // edge_from_v is v -> u
                     int u = edge_from_v.to;
                     int rev_edge_idx = edge_from_v.rev;

                     // --- Safety Check ---
                     // Ensure u and rev_edge_idx are valid before accessing graph[u]
                     if (u < 0 || u >= n || rev_edge_idx < 0 || rev_edge_idx >= graph[u].size()) {
                         // This indicates a potential graph inconsistency. Log or handle.
                         #pragma omp critical (GraphError)
                         {
                            std::cerr << "Warning (GR Parallel): Invalid reverse edge index encountered."
                                      << " From V=" << v << ", To U=" << u << ", RevIdx=" << rev_edge_idx
                                      << ", Graph[u].size=" << (u >= 0 && u < n ? graph[u].size() : -1) << std::endl;
                         }
                         continue;
                     }
                     // --- End Safety Check ---

                     const auto& edge_u_v = graph[u][rev_edge_idx]; // The actual edge u -> v

                     // Check if there's residual capacity FROM u TO v
                     if (edge_u_v.capacity - edge_u_v.flow > 0) {
                         // Check if u has already been visited/updated (height[u] < n)
                         // Use atomic compare-and-swap (CAS) to ensure only one thread
                         // updates height[u] and adds it to the next frontier *for this level*.

                         int expected_height = n;         // Expect 'infinity' (n)
                         int desired_height = current_level + 1; // New height for this level
                         bool updated_by_this_thread = false;

                         // --- Atomic Operation: Check and Set Height ---
                         // We need to atomically check if height[u] == n and, if so, set it to desired_height.
                         // OpenMP 5.0+ provides 'omp atomic compare capture'.
                         // If using older versions, a critical section or std::atomic might be needed.

                         // Using OpenMP atomic compare capture (preferred if available)
                         #if _OPENMP >= 201811 // Check for OpenMP 5.0 or later
                         #pragma omp atomic compare capture // seq_cst memory order by default
                         { if (height[u] == expected_height) { height[u] = desired_height; updated_by_this_thread = true; } }
                         #else
                         // Fallback using a critical section (less performant)
                         #pragma omp critical (UpdateHeightBFS)
                         {
                              if (height[u] == expected_height) { // Double-check inside critical
                                  height[u] = desired_height;
                                  updated_by_this_thread = true;
                              }
                         }
                         #endif
                         // --- End Atomic Operation ---

                         // If this thread performed the update, add 'u' to its local list
                         if (updated_by_this_thread) {
                             local_next_frontiers[tid].push_back(u);
                         }
                     }
                 } // End loop over neighbors of v
            } // End parallel for loop over frontier (implicit barrier here)
        } // End parallel region

        // 4. Combine local frontiers into the global next frontier (Sequential Part)
        // This merge step is typically fast relative to the parallel exploration.
        size_t total_next_size = 0;
        for(int i=0; i<num_threads; ++i) {
            total_next_size += local_next_frontiers[i].size();
        }
        next_frontier.reserve(total_next_size); // Pre-allocate exact size

        for (int i = 0; i < num_threads; ++i) {
             if (!local_next_frontiers[i].empty()) {
                next_frontier.insert(next_frontier.end(),
                                      std::make_move_iterator(local_next_frontiers[i].begin()),
                                      std::make_move_iterator(local_next_frontiers[i].end()));
                local_next_frontiers[i].clear(); // Clear local list after moving
             }
        }


        // Prepare for the next level
        frontier.swap(next_frontier); // Use the newly populated list as the next frontier
        // next_frontier is now empty (or holds old frontier data) and ready for next iteration

        if (!frontier.empty()) {
            current_level++;
        }
    } // End while loop (BFS finished)

    // 5. Finalization
    // Ensure source height is n (unless source == sink, which is handled by init)
    if (n > 0 && source != sink) {
        height[source] = n;
    }

    // --- Stop Timing (Optional) ---
    // #ifdef _OPENMP
    // double end_time = omp_get_wtime();
    // std::cout << "Parallel Global Relabel finished in " << (end_time - start_time)
    //           << " seconds. Levels processed: " << current_level << std::endl;
    // #endif
    // --- End Timing ---
}

// Debug print state (unchanged)
void PushRelabelParallel::printState(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                     const std::vector<int>& excess, const std::vector<int>& height,
                                     const std::queue<int>& active_vertices, int n) {
    std::cout << "--- Current State ---" << std::endl;
    std::cout << "  Heights: ";
    for (int i = 0; i < n; i++) {
        std::cout << (height[i] >= n ? "INF " : std::to_string(height[i]) + " ");
    }
    std::cout << std::endl;
    if (!excess.empty()) {
        std::cout << "  Excess:  ";
        for (int i = 0; i < n; i++) std::cout << excess[i] << " ";
        std::cout << std::endl;
    }
    std::cout << "  Active Queue (" << active_vertices.size() << "): ";
    std::queue<int> tmp = active_vertices;
    while (!tmp.empty()) { std::cout << tmp.front() << " "; tmp.pop(); }
    std::cout << std::endl;
    std::cout << "---------------------" << std::endl;
}