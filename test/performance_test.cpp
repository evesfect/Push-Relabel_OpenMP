#include "../src/flow_network.h"
#include "../src/generator.h"
#include "../src/test.h"
#include "../src/sequential.h"
#include "../src/parallel.h"
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>

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
    
    // Benchmark sequential algorithm
    std::cout << "\nSequential Push-Relabel:" << std::endl;
    double seq_total_time = 0.0;
    int seq_flow = 0;
    
    for (int i = 0; i < runs; ++i) {
        // Create a copy of the network for each run
        FlowNetwork network_copy = network;
        
        Timer timer;
        seq_flow = PushRelabelSequential::maxFlow(network_copy, source, sink);
        double elapsed = timer.elapsed();
        seq_total_time += elapsed;
        
        std::cout << "  Run " << (i+1) << ": " << elapsed << " ms, flow = " << seq_flow << std::endl;
        
        // Verify correctness
        bool correct = (seq_flow == expected_flow);
        if (!correct) {
            std::cout << "  ERROR: Sequential flow doesn't match expected flow" << std::endl;
        }
    }
    
    double seq_avg_time = seq_total_time / runs;
    std::cout << "  Average: " << seq_avg_time << " ms" << std::endl;
    
    // Benchmark parallel algorithm with different thread counts
    int max_threads = 8; // Adjust based on your system
    
    #ifdef _OPENMP
    std::cout << "\nParallel Push-Relabel:" << std::endl;
    
    for (int num_threads = 2; num_threads <= max_threads; num_threads *= 2) {
        double par_total_time = 0.0;
        int par_flow = 0;
        
        std::cout << "\n  Using " << num_threads << " threads:" << std::endl;
        
        for (int i = 0; i < runs; ++i) {
            // Create a copy of the network for each run
            FlowNetwork network_copy = network;
            
            Timer timer;
            par_flow = PushRelabelParallel::maxFlow(network_copy, source, sink, num_threads);
            double elapsed = timer.elapsed();
            par_total_time += elapsed;
            
            std::cout << "    Run " << (i+1) << ": " << elapsed << " ms, flow = " << par_flow << std::endl;
            
            // Verify correctness
            bool correct = (par_flow == expected_flow);
            if (!correct) {
                std::cout << "    ERROR: Parallel flow doesn't match expected flow" << std::endl;
            }
        }
        
        double par_avg_time = par_total_time / runs;
        double speedup = seq_avg_time / par_avg_time;
        
        std::cout << "    Average: " << par_avg_time << " ms" << std::endl;
        std::cout << "    Speedup: " << std::fixed << std::setprecision(2) << speedup << "x" << std::endl;
    }
    #else
    std::cout << "\nOpenMP not enabled, skipping parallel benchmark" << std::endl;
    #endif
}

int main() {
    try {
        runBenchmark(50, 0.5, 1, 300, 3);
        runBenchmark(100, 0.3, 1, 200, 3);       
        runBenchmark(1000, 0.1, 1, 100, 3);
        
        std::cout << "\nAll benchmarks completed." << std::endl;
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}