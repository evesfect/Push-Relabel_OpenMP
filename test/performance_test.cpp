#include "../src/flow_network.h"
#include "../src/generator.h"
#include "../src/test.h"
#include "../src/sequential.h"
#include "../src/parallel.h"
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <chrono>

// Test global relabeling in isolation
void testGlobalRelabelingPerformance(int numVertices, double density) {
    std::cout << "\n=== Testing Global Relabeling Performance ===" << std::endl;
    std::cout << "Graph size: " << numVertices << " vertices, density: " << density << std::endl;

    // Generate a random graph
    FlowNetworkGenerator generator(42);  // Fixed seed for reproducibility
    FlowNetwork network = generator.generateRandom(numVertices, density, 1, 100);
    
    int source = 0;
    int sink = numVertices - 1;

    // Make copies for both implementations
    auto& graph = const_cast<std::vector<std::vector<FlowNetwork::Edge>>&>(network.getGraph());
    std::vector<int> height_seq(numVertices);
    std::vector<int> height_par(numVertices);

    // Test sequential global relabeling
    {
        Timer timer("Sequential Global Relabeling");
        PushRelabelSequential::globalRelabel(graph, height_seq, source, sink, numVertices);
    }

    // Test parallel global relabeling
    {
        Timer timer("Parallel Global Relabeling");
        PushRelabelParallel::globalRelabel(graph, height_par, source, sink, numVertices);
    }

    // Verify both implementations produce the same heights
    bool heights_match = true;
    for (int i = 0; i < numVertices && heights_match; ++i) {
        if (height_seq[i] != height_par[i]) {
            heights_match = false;
            std::cout << "Height mismatch at vertex " << i 
                      << ": seq=" << height_seq[i] 
                      << ", par=" << height_par[i] << std::endl;
        }
    }
    
    std::cout << "Heights match: " << (heights_match ? "YES" : "NO") << std::endl;
}

// Test full max-flow algorithm
void testMaxFlowPerformance(int numVertices, double density, int numThreads) {
    std::cout << "\n=== Testing Max Flow Performance ===" << std::endl;
    std::cout << "Graph size: " << numVertices << " vertices, density: " << density 
              << ", threads: " << numThreads << std::endl;

    // Generate a random graph
    FlowNetworkGenerator generator(42);  // Fixed seed for reproducibility
    FlowNetwork network = generator.generateRandom(numVertices, density, 1, 100);
    
    int source = 0;
    int sink = numVertices - 1;

    // Test sequential implementation
    FlowNetwork network_seq = network;
    int sequential_flow;
    {
        Timer timer("Sequential Max Flow");
        sequential_flow = PushRelabelSequential::maxFlow(network_seq, source, sink);
    }

    // Test parallel implementation
    FlowNetwork network_par = network;
    int parallel_flow;
    {
        Timer timer("Parallel Max Flow");
        parallel_flow = PushRelabelParallel::maxFlow(network_par, source, sink, numThreads);
    }

    std::cout << "Sequential flow: " << sequential_flow << std::endl;
    std::cout << "Parallel flow: " << parallel_flow << std::endl;
    std::cout << "Flows match: " << (sequential_flow == parallel_flow ? "YES" : "NO") << std::endl;
}

int main(int argc, char* argv[]) {
    // Default values
    int numVertices = 20000;
    double density = 0.3;
    int numThreads = 4;
    
    // Parse command line arguments if provided
    if (argc > 1) numVertices = std::stoi(argv[1]);
    if (argc > 2) density = std::stod(argv[2]);
    if (argc > 3) numThreads = std::stoi(argv[3]);
    
    // Test small, medium and large graphs
    std::vector<int> graphSizes = {1000, 5000, 10000};
    std::vector<double> densities = {0.01, 0.1, 0.3, 0.5};
    
    // Test isolated global relabeling performance
    
    for (int size : graphSizes) {
        for (double dens : densities) {
            testGlobalRelabelingPerformance(size, dens);
        }
    
    }
    // Test full max-flow algorithm performance
    //testMaxFlowPerformance(numVertices, density, numThreads);
    
    return 0;
}