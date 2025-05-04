// parallel.h
#ifndef PARALLEL_H
#define PARALLEL_H

#include "flow_network.h"
#include <vector>
#include <iostream>
#include <queue>
#include <omp.h> // Keep OpenMP include for future parallelization
#include <limits>
#include <atomic> // Include for atomic operations

class PushRelabelParallel {
public:
    static const int MAX_ITERATIONS = 100000000; // Consider adjusting if needed
    static const int INF_HEIGHT = std::numeric_limits<int>::max(); // For global relabeling

    static int maxFlow(FlowNetwork& network, int source, int sink, int num_threads = 0);
    static int maxFlow_ActiveParallel(FlowNetwork& network, int source, int sink, int num_threads = 0);
    static void globalRelabel(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                             std::vector<int>& height, int source, int sink, int n);
    static void globalRelabel_timed(std::vector<std::vector<FlowNetwork::Edge>>& graph, // Function with detailed timings
                                    std::vector<int>& height, int source, int sink, int n);

private:
    // Initialization
    static void initialize(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                          std::vector<int>& excess, std::vector<int>& height,
                          std::queue<int>& active_vertices, std::vector<bool>& in_queue, // Changed to vector<bool>
                          int source, int sink, int n);

    // Push operation
    static bool push(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                    std::vector<int>& excess, std::vector<int>& height, // Added height
                    std::queue<int>& active_vertices, std::vector<bool>& in_queue,
                    int u, int v_idx, int source, int sink); // Added height, source, sink

    // Relabel operation
    static bool relabel(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                       std::vector<int>& height, int u, int n); // Added n

    // Discharge operation
    static int discharge(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                        std::vector<int>& excess, std::vector<int>& height,
                        std::queue<int>& active_vertices, std::vector<bool>& in_queue,
                        int u, int source, int sink, int n); // Added n, returns relabel count

    // New discharge variant using atomic operations (called by parallel loop)
    // It collects newly active nodes locally instead of modifying the global queue directly.
    static int discharge_Atomic(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                std::vector<std::atomic<int>>& atomic_excess, // Use atomic wrapper
                                std::vector<int>& height,
                                std::vector<int>& local_newly_active, // Thread-local list
                                int u, int source, int sink, int n);

    // Debug helper
    static void printState(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                          const std::vector<int>& excess, const std::vector<int>& height,
                          const std::queue<int>& active_vertices, int n);
     // Overload for atomic excess if needed for debugging
    static void printState(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                          const std::vector<std::atomic<int>>& atomic_excess, const std::vector<int>& height,
                          const std::queue<int>& active_vertices, int n);
};

#endif // PARALLEL_H