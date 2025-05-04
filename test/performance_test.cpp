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
    FlowNetwork network_seq = network; // Copy for sequential
    int sequential_flow;
    {
        Timer timer("Sequential Max Flow");
        sequential_flow = PushRelabelSequential::maxFlow(network_seq, source, sink);
    }

    // Test original parallel implementation
    FlowNetwork network_par_orig = network; // Copy for original parallel
    int parallel_flow_orig;
    {
        Timer timer("Parallel Max Flow (Original)");
        parallel_flow_orig = PushRelabelParallel::maxFlow(network_par_orig, source, sink, numThreads);
    }

    // Test new active-parallel implementation
    FlowNetwork network_par_active = network; // Copy for active-parallel
    int parallel_flow_active;
    {
        Timer timer("Parallel Max Flow (ActiveParallel)");
        parallel_flow_active = PushRelabelParallel::maxFlow_ActiveParallel(network_par_active, source, sink, numThreads);
    }

    std::cout << "--- Results ---" << std::endl;
    std::cout << "Sequential flow:           " << sequential_flow << std::endl;
    std::cout << "Parallel flow (Original):  " << parallel_flow_orig << std::endl;
    std::cout << "Parallel flow (ActivePar): " << parallel_flow_active << std::endl;
    
    std::cout << "\n--- Comparisons ---" << std::endl;
    bool seq_vs_orig = (sequential_flow == parallel_flow_orig);
    bool seq_vs_active = (sequential_flow == parallel_flow_active);
    bool orig_vs_active = (parallel_flow_orig == parallel_flow_active);

    std::cout << "Flows match (Seq vs Orig):    " << (seq_vs_orig ? "YES" : "NO") << std::endl;
    std::cout << "Flows match (Seq vs Active):  " << (seq_vs_active ? "YES" : "NO") << std::endl;
    std::cout << "Flows match (Orig vs Active): " << (orig_vs_active ? "YES" : "NO") << std::endl;

    if (!seq_vs_orig || !seq_vs_active) {
         std::cerr << "WARNING: Flow mismatch detected!" << std::endl;
    }
}

int main(int argc, char* argv[]) {
    // Default values
    int numVertices = 15000;
    double density = 0.3;
    int numThreads = 4;
    
    // Parse command line arguments if provided
    if (argc > 1) numVertices = std::stoi(argv[1]);
    if (argc > 2) density = std::stod(argv[2]);
    if (argc > 3) numThreads = std::stoi(argv[3]);
    
    // Test small, medium and large graphs for global relabel
    std::vector<int> graphSizesGR = {1000, 5000}; // Reduced for quicker testing
    std::vector<double> densitiesGR = {0.1, 0.3};  // Reduced for quicker testing
    /*
    std::cout << "--- Global Relabel Tests ---" << std::endl;
    for (int size : graphSizesGR) {
        for (double dens : densitiesGR) {
            testGlobalRelabelingPerformance(size, dens);
        }
    }*/
    
    std::cout << "\n--- Max Flow Tests ---" << std::endl;
    // Test full max-flow algorithm performance using specified parameters
    testMaxFlowPerformance(numVertices, density, numThreads);
    
    return 0;
}