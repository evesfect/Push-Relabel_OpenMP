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
    
    // Test our sequential implementation
    int sequential_flow = PushRelabelSequential::maxFlow(network, source, sink);
    bool correct = MaxFlowTester::verifyCorrectness(network, source, sink, sequential_flow);
    
    if (correct) {
        std::cout << "Tiny graph test passed!" << std::endl;
    } else {
        std::cout << "Tiny graph test failed!" << std::endl;
    }
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
    
    // Test our sequential implementation
    int sequential_flow = PushRelabelSequential::maxFlow(network, source, sink);
    bool correct = MaxFlowTester::verifyCorrectness(network, source, sink, sequential_flow);
    
    if (correct) {
        std::cout << "Small graph test passed!" << std::endl;
    } else {
        std::cout << "Small graph test failed!" << std::endl;
    }
}

void testRandomGraph() {
    std::cout << "\n=== Testing Random Graph ===" << std::endl;
    
    // Generate a random flow network (smaller than before)
    FlowNetworkGenerator generator(42); // Fixed seed for reproducibility
    FlowNetwork network = generator.generateRandom(10, 0.3, 1, 100);
    
    // Source and sink
    int source = 0;
    int sink = 9;
    
    // Compute expected max flow using Boost
    int expected_flow = MaxFlowTester::computeExpectedFlow(network, source, sink);
    std::cout << "Expected max flow: " << expected_flow << std::endl;
    
    // Test our sequential implementation
    int sequential_flow = PushRelabelSequential::maxFlow(network, source, sink);
    bool correct = MaxFlowTester::verifyCorrectness(network, source, sink, sequential_flow);
    
    if (correct) {
        std::cout << "Random graph test passed!" << std::endl;
    } else {
        std::cout << "Random graph test failed!" << std::endl;
    }
}

int main() {
    try {
        testTinyGraph();
        testSmallGraph();
        testRandomGraph();
        
        std::cout << "\nTesting terminated." << std::endl;
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}