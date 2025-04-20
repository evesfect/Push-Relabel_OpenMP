#ifndef PARALLEL_H
#define PARALLEL_H

#include "flow_network.h"
#include <vector>
#include <queue>

// Only include OpenMP header if OpenMP is enabled
#ifdef _OPENMP
#include <omp.h>
#endif

class PushRelabelParallel {
public:
    // Maximum number of iterations to prevent infinite loops
    static const int MAX_ITERATIONS = 1000000;
    
    // Find maximum flow using parallel push-relabel algorithm
    static int maxFlow(FlowNetwork& network, int source, int sink, int num_threads = 0);

private:
    // Batch size for queue operations to reduce contention
    static const int BATCH_SIZE = 32;
    
    // Threshold for global relabeling (operations per vertex)
    static const int GLOBAL_RELABEL_FREQ = 2;
    
    // Initialize preflow from source
    static void initialize(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                          std::vector<int>& excess, std::vector<int>& height,
                          std::queue<int>& active_vertices, bool* in_queue,
                          int source, int sink, int n);
    
    // Discharge operation with reduced atomics
    static void discharge(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                          std::vector<int>& excess, std::vector<int>& height,
                          std::vector<int>& local_to_add, bool* in_queue,
                          int u, int source, int sink);
    
    // Global relabeling using parallel backward BFS
    static void globalRelabeling(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                std::vector<int>& height, 
                                std::vector<int>& excess,
                                std::queue<int>& active_vertices,
                                bool* in_queue,
                                int source, int sink, int n);
    
    // Dequeue a batch of vertices
    static void dequeueBatch(std::queue<int>& global_queue, 
                            std::vector<int>& local_batch,
                            int max_size);
    
    // Add a batch of vertices to the queue
    static void enqueueBatch(std::queue<int>& global_queue,
                           std::vector<int>& local_batch,
                           bool* in_queue);
};

#endif // PARALLEL_H