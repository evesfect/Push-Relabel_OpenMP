#include "parallel.h"
#include "sequential.h"
#include <limits>
#include <algorithm>
#include <iomanip>
#include <vector>
#include <atomic>
#include <numeric>

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

    int relabel_since_last = 0;
    const int global_freq = n;
    int total_global = 0;

    initialize(graph, excess, height, active_vertices, in_queue, source, sink, n);

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

// Global relabel via backward BFS version 2 (TIMED VERSION) (lots of overhead!)
void PushRelabelParallel::globalRelabel_timed(std::vector<std::vector<FlowNetwork::Edge>>& graph, // Needs non-const graph
                                              std::vector<int>& height, int source, int sink, int n) {
    double start_time_total = 0.0;
    #ifdef _OPENMP
    start_time_total = omp_get_wtime();
    #endif

    double start_time_init = 0.0;
    #ifdef _OPENMP
    start_time_init = omp_get_wtime();
    #endif

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        height[i] = n;
    }
    height[sink] = 0;

    double end_time_init = 0.0;
    #ifdef _OPENMP
    end_time_init = omp_get_wtime();
    #endif

    double start_time_setup = 0.0;
    #ifdef _OPENMP
    start_time_setup = omp_get_wtime();
    #endif

    std::vector<int> frontier;
    frontier.reserve(n/4);
    frontier.push_back(sink);

    std::vector<int> next_frontier;
    next_frontier.reserve(n/4);

    int current_level = 0;
    int num_threads = 1;
    #ifdef _OPENMP
    num_threads = omp_get_max_threads();
    #endif

    std::vector<std::vector<int>> local_frontiers(num_threads);
    for (int i = 0; i < num_threads; ++i) {
        local_frontiers[i].reserve(n/8); // Reserve space but avoid excessive allocation
    }

    double end_time_setup = 0.0;
    #ifdef _OPENMP
    end_time_setup = omp_get_wtime();
    #endif

    double start_time_bfs = 0.0;
    double total_time_merge = 0.0;
    #ifdef _OPENMP
    start_time_bfs = omp_get_wtime();
    #endif

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
        }

        // Merge timing
        double start_time_merge = 0.0;
        #ifdef _OPENMP
        start_time_merge = omp_get_wtime();
        #endif

        // Merge local frontiers
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


    #ifdef _OPENMP
    std::cout << "Global Relabel Timings:" << std::endl
              << "  - Initialization: " << (end_time_init - start_time_init) * 1000 << " ms" << std::endl
              << "  - Setup:         " << (end_time_setup - start_time_setup) * 1000 << " ms" << std::endl
              << "  - BFS Loop:      " << (end_time_bfs - start_time_bfs) * 1000 << " ms" << std::endl
              << "  - Merge Time:  " << total_time_merge * 1000 << " ms" << std::endl
              << "  - Finalization:  " << (end_time_finalize - start_time_finalize) * 1000 << " ms" << std::endl
              << "  - Total Time:    " << (end_time_total - start_time_total) * 1000 << " ms" << std::endl
              << "  - Levels:        " << current_level << std::endl;
    #endif
}

// Global relabel Hybrid
void PushRelabelParallel::globalRelabel(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                        std::vector<int>& height, int source, int sink, int n) {

    const int SEQUENTIAL_THRESHOLD_VERTICES = 1300;

    if (n < SEQUENTIAL_THRESHOLD_VERTICES) {
        PushRelabelSequential::globalRelabel(graph, height, source, sink, n);
    } else {
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i) {
            height[i] = n;
        }
        if (n > 0) { height[sink] = 0; }

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
            }

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

        if (n > 0 && source != sink) {
            height[source] = n;
        }
    }
}

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


int PushRelabelParallel::discharge_Atomic(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                          std::vector<std::atomic<int>>& atomic_excess,
                                          std::vector<int>& height,
                                          std::vector<int>& local_newly_active,
                                          int u, int source, int sink, int n) {
    if (u == source || u == sink) return 0;

    int relabels = 0;
    int current_excess = atomic_excess[u].load(std::memory_order_relaxed); // Load initial excess

    while (current_excess > 0 && height[u] < n) {
        bool did_push = false;
        for (size_t i = 0; i < graph[u].size() && current_excess > 0; ++i) {
            auto& e = graph[u][i];
            int v = e.to;

            // Check if push is possible
            if (e.capacity - e.flow > 0 && height[u] == height[v] + 1) {
                auto& rev = graph[v][e.rev];
                int delta = std::min(current_excess, e.capacity - e.flow);

                if (delta > 0) {
                    int expected_u = current_excess;
                    while (expected_u > 0 && !atomic_excess[u].compare_exchange_weak(expected_u, expected_u - delta)) {
                         delta = std::min(expected_u, e.capacity - e.flow);
                         if (delta <= 0) break; // If excess becomes 0 or capacity constraint changes
                    }

                    if (delta > 0 && expected_u > 0) { // If compare_exchange succeeded and delta 
                        #pragma omp atomic update
                        e.flow += delta;
                        #pragma omp atomic update
                        rev.flow -= delta;

                        atomic_excess[v].fetch_add(delta, std::memory_order_relaxed);

                        did_push = true;
                        current_excess -= delta; // Update local view of excess

                        // If v becomes active and is not source/sink, then add to local list
                        if (v != source && v != sink && atomic_excess[v].load(std::memory_order_relaxed) > 0) {
                             local_newly_active.push_back(v);
                        }
                    } else {
                        // Failed 
                        current_excess = atomic_excess[u].load(std::memory_order_relaxed);
                        did_push = false;
                        break; // re-evaluaete
                    }
                }
            }
             // Reload excess if needed, especially if failed
            if (!did_push) current_excess = atomic_excess[u].load(std::memory_order_relaxed);
        }

        if (!did_push && current_excess > 0) {
            // Relabel if no push but excess remains
            if (relabel(graph, height, u, n)) {
                relabels++;
            } else {
                break;
            }
        }
        current_excess = atomic_excess[u].load(std::memory_order_relaxed);

    }

    return relabels;
}

int PushRelabelParallel::maxFlow_ActiveParallel(FlowNetwork& network, int source, int sink, int num_threads) {
     int n = network.getNumVertices();
    std::cout << "Starting Parallel Push-Relabel (ActiveParallel) with " << n
              << " vertices (source=" << source << ", sink=" << sink << ")" << std::endl;

    #ifdef _OPENMP
    if (num_threads <= 0) {
        num_threads = omp_get_max_threads();
    }
    omp_set_num_threads(num_threads);
    std::cout << "Using OpenMP threads: " << num_threads << std::endl;
    #else
    num_threads = 1;
    std::cout << "Warning: OpenMP not enabled, running sequentially." << std::endl;
    #endif

    if (source < 0 || source >= n || sink < 0 || sink >= n || source == sink) {
        throw std::invalid_argument("Invalid source or sink");
    }

    auto& graph = const_cast<std::vector<std::vector<FlowNetwork::Edge>>&>(network.getGraph());

    std::vector<int> initial_excess(n, 0);
    std::vector<int> height(n, 0);
    std::queue<int> active_vertices_queue;
    std::vector<bool> in_queue(n, false);

    initialize(graph, initial_excess, height, active_vertices_queue, in_queue, source, sink, n);

    std::vector<std::atomic<int>> atomic_excess(n);
    for (int i = 0; i < n; ++i) {
        atomic_excess[i].store(initial_excess[i], std::memory_order_relaxed);
    }
    initial_excess.clear();


    int relabel_since_last = 0;
    const int global_freq = n;
    int total_global = 0;

    globalRelabel(graph, height, source, sink, n);
    total_global++;

    active_vertices_queue = {};
    std::fill(in_queue.begin(), in_queue.end(), false);
    for (int i = 0; i < n; ++i) {
        if (i != source && i != sink && atomic_excess[i].load(std::memory_order_relaxed) > 0 && height[i] < n) {
            active_vertices_queue.push(i);
            in_queue[i] = true;
        }
    }

    int iterations = 0, total_local_relabels = 0;
    std::vector<int> current_active_batch;
    std::vector<std::vector<int>> thread_local_newly_active(num_threads);


    while (!active_vertices_queue.empty() && iterations < MAX_ITERATIONS) {
        if (relabel_since_last >= global_freq) {
            globalRelabel(graph, height, source, sink, n);
            total_global++;
            relabel_since_last = 0;

            active_vertices_queue = {};
            std::fill(in_queue.begin(), in_queue.end(), false);
            for (int i = 0; i < n; ++i) {
                 if (i != source && i != sink && atomic_excess[i].load(std::memory_order_relaxed) > 0 && height[i] < n) {
                    active_vertices_queue.push(i);
                    in_queue[i] = true;
                }
            }
            if (active_vertices_queue.empty()) break;
        }

        // Prepare batch
        current_active_batch.clear();
        while (!active_vertices_queue.empty()) {
            current_active_batch.push_back(active_vertices_queue.front());
            active_vertices_queue.pop();
        }

        if (current_active_batch.empty()) continue; // safety check


        // Parallel discharge
        int batch_relabels = 0;

        #pragma omp parallel num_threads(num_threads) reduction(+:batch_relabels) // Not: Combines values of each thread*
        {
            int tid = omp_get_thread_num();
            thread_local_newly_active[tid].clear();

            #pragma omp for schedule(dynamic, 64)
            for (size_t i = 0; i < current_active_batch.size(); ++i) {
                int u = current_active_batch[i];
                if (atomic_excess[u].load(std::memory_order_relaxed) > 0 && height[u] < n) {
                     batch_relabels += discharge_Atomic(graph, atomic_excess, height,
                                                     thread_local_newly_active[tid],
                                                     u, source, sink, n);
                }
                #pragma omp critical (UpdateInQueue) // Protect access to shared in_queue
                {
                in_queue[u] = false;
                }

            }
        }

        total_local_relabels += batch_relabels;
        relabel_since_last += batch_relabels;

        // Merge thread-local newly active nodes into the global queue
        for (int tid = 0; tid < num_threads; ++tid) {
            for (int v : thread_local_newly_active[tid]) {
                if (!in_queue[v] && atomic_excess[v].load(std::memory_order_relaxed) > 0 && height[v] < n) {
                     active_vertices_queue.push(v);
                    in_queue[v] = true;
                }
            }
        }
        for (int u : current_active_batch) {
            if (!in_queue[u] && atomic_excess[u].load(std::memory_order_relaxed) > 0 && height[u] < n) {
                active_vertices_queue.push(u);
                in_queue[u] = true;
            }
        }

        iterations++;
    }

    if (iterations >= MAX_ITERATIONS) {
        std::cout << "WARNING: Reached maximum iterations. Max flow might be incorrect." << std::endl;
    }

    int max_flow = atomic_excess[sink].load(std::memory_order_relaxed);
    std::cout << "Completed (ActiveParallel) in " << iterations << " iterations. Total local relabels: "
              << total_local_relabels << ", global relabels: " << total_global
              << ". Max flow: " << max_flow << std::endl;
    return max_flow;
}

void PushRelabelParallel::printState(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                     const std::vector<std::atomic<int>>& atomic_excess,
                                     const std::vector<int>& height,
                                     const std::queue<int>& active_vertices, int n) {
    std::cout << "--- Current State (Atomic Excess) ---" << std::endl;
    std::cout << "  Heights: ";
    for (int i = 0; i < n; i++) {
        std::cout << (height[i] >= n ? "INF " : std::to_string(height[i]) + " ");
    }
    std::cout << std::endl;
    if (!atomic_excess.empty()) {
        std::cout << "  Excess:  ";
        for (int i = 0; i < n; i++) std::cout << atomic_excess[i].load() << " ";
        std::cout << std::endl;
    }
    std::cout << "  Active Queue (" << active_vertices.size() << "): ";
    std::queue<int> tmp = active_vertices;
    while (!tmp.empty()) { std::cout << tmp.front() << " "; tmp.pop(); }
    std::cout << std::endl;
    std::cout << "-------------------------------------" << std::endl;
}