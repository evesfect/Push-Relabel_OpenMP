#include "../src/flow_network.h"
#include "../src/generator.h"
#include "../src/test.h"
#include "../src/sequential.h"
#include "../src/parallel.h"
#include <iostream>
#include <cassert>

void testTinyGraph() {
    std::cout << "=== Testing Tiny Graph ===" << std::endl;
    
    // Create a very small flow network manually
    FlowNetwork network(4);
    
    // Add edges (source=0, sink=3)
    network.addEdge(0, 1, 3);
    network.addEdge(0, 2, 2);
    network.addEdge(1, 2, 1);
    network.addEdge(1, 3, 3);
    network.addEdge(2, 3, 2);
    
    // Source and sink
    int source = 0;
    int sink = 3;
    
    // Print the network
    network.print();
    
    // Compute expected max flow using Boost
    int expected_flow = MaxFlowTester::computeExpectedFlow(network, source, sink);
    std::cout << "Expected max flow: " << expected_flow << std::endl;
    
    // Test sequential implementation
    FlowNetwork network_seq = network;
    int sequential_flow = PushRelabelSequential::maxFlow(network_seq, source, sink);
    bool seq_correct = MaxFlowTester::verifyCorrectness(network_seq, source, sink, sequential_flow);
    
    std::cout << "Sequential result: " << (seq_correct ? "CORRECT" : "INCORRECT") 
              << " (flow = " << sequential_flow << ")" << std::endl;
    
    // Test parallel implementation
    #ifdef _OPENMP
    FlowNetwork network_par = network;
    int parallel_flow = PushRelabelParallel::maxFlow(network_par, source, sink, 2);
    bool par_correct = MaxFlowTester::verifyCorrectness(network_par, source, sink, parallel_flow);
    
    std::cout << "Parallel result: " << (par_correct ? "CORRECT" : "INCORRECT") 
              << " (flow = " << parallel_flow << ")" << std::endl;
    #else
    std::cout << "OpenMP not enabled, skipping parallel test" << std::endl;
    #endif
}

void testSmallGraph() {
    std::cout << "\n=== Testing Small Graph ===" << std::endl;
    
    // Create a small flow network manually
    FlowNetwork network(6);
    
    // Add edges
    network.addEdge(0, 1, 16);
    network.addEdge(0, 2, 13);
    network.addEdge(1, 3, 12);
    network.addEdge(2, 1, 4);
    network.addEdge(2, 4, 14);
    network.addEdge(3, 2, 9);
    network.addEdge(3, 5, 20);
    network.addEdge(4, 3, 7);
    network.addEdge(4, 5, 4);
    
    // Source and sink
    int source = 0;
    int sink = 5;
    
    // Print the network
    network.print();
    
    // Compute expected max flow using Boost
    int expected_flow = MaxFlowTester::computeExpectedFlow(network, source, sink);
    std::cout << "Expected max flow: " << expected_flow << std::endl;
    
    // Test sequential implementation
    FlowNetwork network_seq = network;
    int sequential_flow = PushRelabelSequential::maxFlow(network_seq, source, sink);
    bool seq_correct = MaxFlowTester::verifyCorrectness(network_seq, source, sink, sequential_flow);
    
    std::cout << "Sequential result: " << (seq_correct ? "CORRECT" : "INCORRECT") 
              << " (flow = " << sequential_flow << ")" << std::endl;
    
    // Test parallel implementation
    #ifdef _OPENMP
    FlowNetwork network_par = network;
    int parallel_flow = PushRelabelParallel::maxFlow(network_par, source, sink, 2);
    bool par_correct = MaxFlowTester::verifyCorrectness(network_par, source, sink, parallel_flow);
    
    std::cout << "Parallel result: " << (par_correct ? "CORRECT" : "INCORRECT") 
              << " (flow = " << parallel_flow << ")" << std::endl;
    #else
    std::cout << "OpenMP not enabled, skipping parallel test" << std::endl;
    #endif
}

void testRandomGraph() {
    std::cout << "\n=== Testing Random Graph ===" << std::endl;
    
    // Generate a random flow network (smaller for testing)
    FlowNetworkGenerator generator(42); // Fixed seed for reproducibility
    FlowNetwork network = generator.generateRandom(10, 0.3, 1, 100);
    
    // Source and sink
    int source = 0;
    int sink = 9;
    
    // Compute expected max flow using Boost
    int expected_flow = MaxFlowTester::computeExpectedFlow(network, source, sink);
    std::cout << "Expected max flow: " << expected_flow << std::endl;
    
    // Test sequential implementation
    FlowNetwork network_seq = network;
    int sequential_flow = PushRelabelSequential::maxFlow(network_seq, source, sink);
    bool seq_correct = MaxFlowTester::verifyCorrectness(network_seq, source, sink, sequential_flow);
    
    std::cout << "Sequential result: " << (seq_correct ? "CORRECT" : "INCORRECT") 
              << " (flow = " << sequential_flow << ")" << std::endl;
    
    // Test parallel implementation
    #ifdef _OPENMP
    FlowNetwork network_par = network;
    int parallel_flow = PushRelabelParallel::maxFlow(network_par, source, sink, 2);
    bool par_correct = MaxFlowTester::verifyCorrectness(network_par, source, sink, parallel_flow);
    
    std::cout << "Parallel result: " << (par_correct ? "CORRECT" : "INCORRECT") 
              << " (flow = " << parallel_flow << ")" << std::endl;
    #else
    std::cout << "OpenMP not enabled, skipping parallel test" << std::endl;
    #endif
}

int main() {
    try {
        // Start with the smallest test
        testTinyGraph();
        
        // Then test slightly larger graphs
        testSmallGraph();
        testRandomGraph();
        
        std::cout << "\nAll tests completed!" << std::endl;
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}