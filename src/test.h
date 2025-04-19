#ifndef TEST_H
#define TEST_H

#include "flow_network.h"
#include <string>
#include <vector>
#include <chrono>

// Forward declarations
class PushRelabelSequential; // Will be implemented later
class PushRelabelParallel;   // Will be implemented later

class MaxFlowTester {
public:
    // Test correctness with Boost's push-relabel algorithm
    static bool verifyCorrectness(const FlowNetwork& network, int source, int sink, int computed_flow);
    
    // Use Boost to compute the expected max flow
    static int computeExpectedFlow(const FlowNetwork& network, int source, int sink);
    
    // Compare sequential and parallel implementations
    static bool compareImplementations(const FlowNetwork& network, int source, int sink);

    // Run benchmarks and print results
    static void runBenchmark(const std::vector<FlowNetwork>& networks, 
                            int source, int sink,
                            int num_runs = 5);
};

// A simple timer class
class Timer {
public:
    Timer() : start_time(std::chrono::high_resolution_clock::now()) {}
    
    void reset() {
        start_time = std::chrono::high_resolution_clock::now();
    }
    
    // Returns elapsed time in milliseconds
    double elapsed() const {
        auto now = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double, std::milli>(now - start_time).count();
    }
    
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
};

#endif // TEST_H