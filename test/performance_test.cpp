#include "../src/flow_network.h"
#include "../src/generator.h"
#include "../src/test.h"
#include "../src/sequential.h"
#include "../src/parallel.h"
#include <iostream>
#include <vector>
#include <string>

void runBenchmark(int size, double edge_prob, int min_cap, int max_cap, int runs) {
    std::cout << "=== Benchmark: Network size " << size 
              << ", edge probability " << edge_prob << " ===" << std::endl;
    
    // Generate a test network
    FlowNetworkGenerator generator;
    FlowNetwork network = generator.generateRandom(size, edge_prob, min_cap, max_cap);
    
    int source = 0;
    int sink = size - 1;
    
    // Verify the expected flow
    int expected_flow = MaxFlowTester::computeExpectedFlow(network, source, sink);
    std::cout << "Expected max flow: " << expected_flow << std::endl;
    
    // Time the sequential algorithm
    Timer timer;
    double total_time = 0.0;
    
    std::cout << "Sequential Push-Relabel:" << std::endl;
    for (int i = 0; i < runs; ++i) {
        network.resetFlow();
        timer.reset();
        
        // Sequential implementation will be added later
        // int flow = PushRelabelSequential::maxFlow(network, source, sink);
        
        double elapsed = timer.elapsed();
        total_time += elapsed;
        
        // std::cout << "  Run " << (i+1) << ": " << elapsed << " ms" << std::endl;
    }
    
    std::cout << "  Average: " << (total_time / runs) << " ms" << std::endl;
    
    // Time the parallel algorithm (will be added later)
    std::cout << "Parallel Push-Relabel:" << std::endl;
    std::cout << "  Not yet implemented" << std::endl;
}

int main() {
    try {
        // Run benchmarks for different sized networks
        runBenchmark(100, 0.1, 1, 100, 5);
        runBenchmark(500, 0.05, 1, 100, 5);
        runBenchmark(1000, 0.01, 1, 100, 3);
        
        std::cout << "\nAll benchmarks completed." << std::endl;
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}