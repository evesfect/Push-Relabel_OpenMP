#ifndef TEST_H
#define TEST_H

#include "flow_network.h"
#include <string>
#include <vector>
#include <chrono>
#include <iostream>

// Forward declarations
class PushRelabelSequential; // Will be implemented later
class PushRelabelParallel;   // Will be implemented later

class MaxFlowTester {
public:
    // Compare results with Boost
    static bool verifyCorrectness(const FlowNetwork& network, int source, int sink, int computed_flow);
    
    // Compute max flow using Boost
    static int computeExpectedFlow(const FlowNetwork& network, int source, int sink);
    
    // Compare sequential and parallel implementations
    static bool compareImplementations(const FlowNetwork& network, int source, int sink);

    // Run benchmarks, print results
    static void runBenchmark(const std::vector<FlowNetwork>& networks, 
                            int source, int sink,
                            int num_runs = 5);
};

// A simple timer class
class Timer {
private:
    std::chrono::high_resolution_clock::time_point startTime;
    std::string name;
    bool print_on_destruction;

public:
    // Constructor to start timing and optionally print on destruction
    Timer(const std::string& timer_name = "", bool print = true)
        : name(timer_name), print_on_destruction(print) {
        startTime = std::chrono::high_resolution_clock::now();
    }

    // Destructor to stop timing and print if requested
    ~Timer() {
        if (print_on_destruction) {
            auto endTime = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
            if (!name.empty()) {
                std::cout << name << ": " << duration << " ms" << std::endl;
            } else {
                std::cout << "Timer: " << duration << " ms" << std::endl;
            }
        }
    }

    void reset() {
        startTime = std::chrono::high_resolution_clock::now();
    }

    // Returns elapsed time in milliseconds without stopping/printing
    double elapsed() const {
        auto now = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double, std::milli>(now - startTime).count();
    }
};

#endif