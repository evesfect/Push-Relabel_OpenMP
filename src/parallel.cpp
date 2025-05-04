#include "parallel.h"
#include "sequential.h"
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

// Global relabel via backward BFS (TIMED VERSION)
void PushRelabelParallel::globalRelabel_timed(std::vector<std::vector<FlowNetwork::Edge>>& graph, // Needs non-const graph
                                              std::vector<int>& height, int source, int sink, int n) {
    // Start overall timing
    double start_time_total = 0.0;
    #ifdef _OPENMP
    start_time_total = omp_get_wtime();
    #endif

    // 1. Initialization timing
    double start_time_init = 0.0;
    #ifdef _OPENMP
    start_time_init = omp_get_wtime();
    #endif

    // Initialize heights array
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        height[i] = n;
    }
    height[sink] = 0;

    double end_time_init = 0.0;
    #ifdef _OPENMP
    end_time_init = omp_get_wtime();
    #endif

    // 2. BFS Setup timing
    double start_time_setup = 0.0;
    #ifdef _OPENMP
    start_time_setup = omp_get_wtime();
    #endif

    // Use vector for frontier (better cache locality than queue)
    std::vector<int> frontier;
    frontier.reserve(n/4); // Reserve reasonable space but not too much
    frontier.push_back(sink);

    std::vector<int> next_frontier;
    next_frontier.reserve(n/4);

    int current_level = 0;
    int num_threads = 1;
    #ifdef _OPENMP
    num_threads = omp_get_max_threads();
    #endif

    // Thread-local frontiers - preallocate once and reuse
    std::vector<std::vector<int>> local_frontiers(num_threads);
    for (int i = 0; i < num_threads; ++i) {
        local_frontiers[i].reserve(n/8); // Reserve space but avoid excessive allocation
    }

    double end_time_setup = 0.0;
    #ifdef _OPENMP
    end_time_setup = omp_get_wtime();
    #endif

    // 3. BFS Loop timing
    double start_time_bfs = 0.0;
    double total_time_merge = 0.0;
    #ifdef _OPENMP
    start_time_bfs = omp_get_wtime();
    #endif

    // Main BFS loop
    while (!frontier.empty()) {
        // Clear all local frontiers at once
        for (auto& local_frontier : local_frontiers) {
            local_frontier.clear();
        }
        next_frontier.clear();

        // Process current level in parallel
        #pragma omp parallel
        {
            int tid = 0;
            #ifdef _OPENMP
            tid = omp_get_thread_num();
            #endif
            auto& my_local_frontier = local_frontiers[tid];

            // Parallel processing of frontier nodes
            #pragma omp for schedule(dynamic, 64)
            for (size_t i = 0; i < frontier.size(); ++i) {
                int v = frontier[i];

                // For each node v, check incoming edges in residual graph
                for (const auto& edge_from_v : graph[v]) {
                    int u = edge_from_v.to;
                    int rev_edge_idx = edge_from_v.rev;

                    // Safety check
                    if (u < 0 || u >= n || rev_edge_idx < 0 || 
                        rev_edge_idx >= static_cast<int>(graph[u].size())) {
                        continue;
                    }

                    const auto& edge_u_v = graph[u][rev_edge_idx];

                    // Check if there's residual capacity from u to v
                    if (edge_u_v.capacity - edge_u_v.flow > 0) {
                        // Try to update height of u using a critical section
                        bool updated = false;
                        // Fast check outside critical section to reduce contention
                        if (height[u] == n) { 
                            #pragma omp critical (UpdateHeightBFS)
                            {
                                // Double-check inside critical section
                                if (height[u] == n) { 
                                    height[u] = current_level + 1;
                                    updated = true;
                                }
                            }
                        }
                        if (updated) {
                            my_local_frontier.push_back(u);
                        }
                    }
                }
            }
        } // End parallel region

        // Merge timing
        double start_time_merge = 0.0;
        #ifdef _OPENMP
        start_time_merge = omp_get_wtime();
        #endif

        // Efficient merge of local frontiers
        size_t total_size = 0;
        for (const auto& local_frontier : local_frontiers) {
            total_size += local_frontier.size();
        }
        next_frontier.reserve(total_size);

        for (const auto& local_frontier : local_frontiers) {
            next_frontier.insert(next_frontier.end(), local_frontier.begin(), local_frontier.end());
        }

        double end_time_merge = 0.0;
        #ifdef _OPENMP
        end_time_merge = omp_get_wtime();
        total_time_merge += (end_time_merge - start_time_merge);
        #endif

        // Prepare for next level
        frontier.swap(next_frontier);
        current_level++;
    }

    double end_time_bfs = 0.0;
    #ifdef _OPENMP
    end_time_bfs = omp_get_wtime();
    #endif

    // 4. Finalization timing
    double start_time_finalize = 0.0;
    #ifdef _OPENMP
    start_time_finalize = omp_get_wtime();
    #endif

    // Set source height to n
    if (source != sink) {
        height[source] = n;
    }

    double end_time_finalize = 0.0;
    double end_time_total = 0.0;
    #ifdef _OPENMP
    end_time_finalize = omp_get_wtime();
    end_time_total = omp_get_wtime();
    #endif

    // Print timing results
    #ifdef _OPENMP
    std::cout << "Global Relabel Timings:" << std::endl
              << "  - Initialization: " << (end_time_init - start_time_init) * 1000 << " ms" << std::endl
              << "  - Setup:         " << (end_time_setup - start_time_setup) * 1000 << " ms" << std::endl
              << "  - BFS Loop:      " << (end_time_bfs - start_time_bfs) * 1000 << " ms" << std::endl
              << "    - Merge Time:  " << total_time_merge * 1000 << " ms" << std::endl
              << "  - Finalization:  " << (end_time_finalize - start_time_finalize) * 1000 << " ms" << std::endl
              << "  - Total Time:    " << (end_time_total - start_time_total) * 1000 << " ms" << std::endl
              << "  - Levels:        " << current_level << std::endl;
    #endif
}

// Global relabel via backward BFS (Hybrid: Sequential for small, Parallel for large)
void PushRelabelParallel::globalRelabel(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                        std::vector<int>& height, int source, int sink, int n) {

    // Define the threshold for switching between sequential and parallel
    const int SEQUENTIAL_THRESHOLD_VERTICES = 1300;

    if (n < SEQUENTIAL_THRESHOLD_VERTICES) {
        // Use the sequential version for small graphs
        PushRelabelSequential::globalRelabel(graph, height, source, sink, n);
        // Optional: Add a print statement here for debugging to see when sequential is used
        // std::cout << "(Using sequential global relabel for n=" << n << ")" << std::endl;
    } else {
        // Use the parallel version for large graphs
        // (Parallel implementation from the previous step goes here)
        // 1. Initialization
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i) {
            height[i] = n;
        }
        if (n > 0) { height[sink] = 0; }

        // 2. BFS Setup
        std::vector<int> frontier;
        if (n > 0) {
            frontier.reserve(n / 4);
            frontier.push_back(sink);
        }
        std::vector<int> next_frontier;
        next_frontier.reserve(n / 4);
        int current_level = 0;
        int num_threads = 1;
        #ifdef _OPENMP
        num_threads = omp_get_max_threads();
        #endif
        std::vector<std::vector<int>> local_frontiers(num_threads);
        for (int i = 0; i < num_threads; ++i) {
            local_frontiers[i].reserve(n / 8);
        }

        // 3. BFS Loop
        while (!frontier.empty()) {
            for (auto& local_frontier : local_frontiers) {
                local_frontier.clear();
            }
            next_frontier.clear();

            #pragma omp parallel
            {
                int tid = 0;
                #ifdef _OPENMP
                tid = omp_get_thread_num();
                #endif
                auto& my_local_frontier = local_frontiers[tid];

                #pragma omp for schedule(dynamic, 64) nowait
                for (size_t i = 0; i < frontier.size(); ++i) {
                    int v = frontier[i];
                    for (const auto& edge_from_v : graph[v]) {
                        int u = edge_from_v.to;
                        int rev_edge_idx = edge_from_v.rev;
                        if (u < 0 || u >= n || rev_edge_idx < 0 || rev_edge_idx >= static_cast<int>(graph[u].size())) {
                            continue;
                        }
                        const auto& edge_u_v = graph[u][rev_edge_idx];
                        if (edge_u_v.capacity - edge_u_v.flow > 0) {
                            bool updated = false;
                            if (height[u] == n) {
                                #pragma omp critical (UpdateHeightBFS)
                                {
                                    if (height[u] == n) {
                                        height[u] = current_level + 1;
                                        updated = true;
                                    }
                                }
                            }
                            if (updated) {
                                my_local_frontier.push_back(u);
                            }
                        }
                    }
                }
            } // End parallel region

            // Merge local frontiers
            size_t total_size = 0;
            for (const auto& local_frontier : local_frontiers) {
                total_size += local_frontier.size();
            }
            if (total_size > next_frontier.capacity()) {
                next_frontier.reserve(total_size);
            }
            for (const auto& local_frontier : local_frontiers) {
                next_frontier.insert(next_frontier.end(), local_frontier.begin(), local_frontier.end());
            }

            frontier.swap(next_frontier);
            if (!frontier.empty()) {
                current_level++;
            }
        }

        // 4. Finalization
        if (n > 0 && source != sink) {
            height[source] = n;
        }
    }
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