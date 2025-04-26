#ifndef PARALLEL_H
#define PARALLEL_H

#include "flow_network.h"
#include <vector>
#include <iostream>
#include <omp.h>

class PushRelabelParallel {
public:
    static const int MAX_ITERATIONS = 100000000;
    
    static int maxFlow(FlowNetwork& network, int source, int sink, int num_threads = 0);
    
private:
    // Initialization
    static void initialize(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                          std::vector<int>& excess, std::vector<int>& height,
                          std::vector<int>& active_vertices, std::vector<bool>& is_active,
                          int& active_count, int source, int sink, int n);
    
    // Phase 1: Push operations - store changes in temporary arrays
    static void pushPhase(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                         const std::vector<int>& excess, const std::vector<int>& height,
                         std::vector<int>& excess_changes, const std::vector<int>& active_vertices,
                         std::vector<bool>& still_active, int active_count, int source, int sink);
    
    // Phase 2: Compute new labels for vertices that are still active
    static void labelComputationPhase(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                     const std::vector<int>& excess, const std::vector<int>& height,
                                     std::vector<int>& new_heights, const std::vector<bool>& still_active);
    
    // Phase 3: Apply the new labels
    static void labelApplicationPhase(std::vector<int>& height, const std::vector<int>& new_heights,
                                     const std::vector<bool>& still_active);
    
    // Phase 4: Apply excess changes and update active vertices
    static void excessUpdatePhase(std::vector<int>& excess, const std::vector<int>& excess_changes,
                                 std::vector<int>& next_active_vertices, std::vector<bool>& is_active,
                                 int& active_count, int source, int sink, int n);
    
    // Helper functions
    static bool push(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                    const std::vector<int>& excess, std::vector<int>& excess_changes,
                    int u, int v_idx);
    
    static bool computeNewLabel(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                              const std::vector<int>& height, int u, int& new_height);
    
    static void printState(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                          const std::vector<int>& excess, const std::vector<int>& height,
                          const std::vector<int>& active_vertices, int active_count, int n);
};

#endif