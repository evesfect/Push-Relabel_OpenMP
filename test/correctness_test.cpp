#include "../src/flow_network.h"
#include "../src/generator.h"
#include "../src/test.h"
#include "../src/sequential.h"
#include "../src/parallel.h"
#include <iostream>
#include <cassert>

void testSmallGraph() {
    std::cout << "=== Testing Small Graph ===" << std::endl;
    
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
    
    // In the future: Test our implementations
    // int sequential_flow = PushRelabelSequential::maxFlow(network, source, sink);
    // MaxFlowTester::verifyCorrectness(network, source, sink, sequential_flow);
}

void testRandomGraph() {
    std::cout << "\n=== Testing Random Graph ===" << std::endl;
    
    // Generate a random flow network
    FlowNetworkGenerator generator(42); // Fixed seed for reproducibility
    FlowNetwork network = generator.generateRandom(20, 0.3, 1, 100);
    
    // Source and sink
    int source = 0;
    int sink = 19;
    
    // Compute expected max flow using Boost
    int expected_flow = MaxFlowTester::computeExpectedFlow(network, source, sink);
    std::cout << "Expected max flow: " << expected_flow << std::endl;
    
    // In the future: Test our implementations
    // int sequential_flow = PushRelabelSequential::maxFlow(network, source, sink);
    // MaxFlowTester::verifyCorrectness(network, source, sink, sequential_flow);
}

void testLayeredGraph() {
    std::cout << "\n=== Testing Layered Graph ===" << std::endl;
    
    // Generate a layered flow network
    FlowNetworkGenerator generator(42); // Fixed seed for reproducibility
    FlowNetwork network = generator.generateLayered(4, 5, 0.5, 1, 100);
    
    // Source and sink
    int source = 0;
    int sink = 19;
    
    // Compute expected max flow using Boost
    int expected_flow = MaxFlowTester::computeExpectedFlow(network, source, sink);
    std::cout << "Expected max flow: " << expected_flow << std::endl;
    
    // In the future: Test our implementations
    // int sequential_flow = PushRelabelSequential::maxFlow(network, source, sink);
    // MaxFlowTester::verifyCorrectness(network, source, sink, sequential_flow);
}

int main() {
    try {
        testSmallGraph();
        testRandomGraph();
        testLayeredGraph();
        
        std::cout << "\nAll tests completed successfully!" << std::endl;
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}