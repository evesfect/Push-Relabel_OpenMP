#include <iostream>
#include <string>
#include <vector>
#include "flow_network.h"
#include "generator.h"
#include "test.h"
#include "sequential.h"
#include "parallel.h"

void printUsage(const std::string& program_name) {
    std::cout << "Usage: " << program_name << " [options]" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -g/--generate <vertices> <edge_prob> <min_cap> <max_cap> <filename>" << std::endl;
    std::cout << "    Generate a random flow network and save it to file" << std::endl;
    std::cout << "  -t/--test <filename> <source> <sink>" << std::endl;
    std::cout << "    Test the max flow algorithm on the given network" << std::endl;
    std::cout << "  -b/--benchmark <filename> <source> <sink> [runs=5]" << std::endl;
    std::cout << "    Benchmark sequential vs parallel implementation" << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printUsage(argv[0]);
        return 1;
    }

    std::string command = argv[1];

    try {
        if (command == "-g" || command == "--generate") {
            if (argc < 7) {
                std::cout << "Error: Not enough arguments for generate command" << std::endl;
                printUsage(argv[0]);
                return 1;
            }

            int vertices = std::stoi(argv[2]);
            double edge_prob = std::stod(argv[3]);
            int min_cap = std::stoi(argv[4]);
            int max_cap = std::stoi(argv[5]);
            std::string filename = argv[6];

            FlowNetworkGenerator generator;
            FlowNetwork network = generator.generateRandom(vertices, edge_prob, min_cap, max_cap);
            
            if (network.saveToFile(filename)) {
                std::cout << "Network saved to " << filename << std::endl;
            } else {
                std::cout << "Error saving network to file" << std::endl;
                return 1;
            }
        }
        else if (command == "-t" || command == "--test") {
            if (argc < 5) {
                std::cout << "Error: Not enough arguments for test command" << std::endl;
                printUsage(argv[0]);
                return 1;
            }

            std::string filename = argv[2];
            int source = std::stoi(argv[3]);
            int sink = std::stoi(argv[4]);

            FlowNetwork network = FlowNetwork::loadFromFile(filename);
            
            // Calculate expected flow using Boost
            int expected_flow = MaxFlowTester::computeExpectedFlow(network, source, sink);
            std::cout << "Expected max flow: " << expected_flow << std::endl;
            
        }
        else if (command == "-b" || command == "--benchmark") {
            if (argc < 5) {
                std::cout << "Error: Not enough arguments for benchmark command" << std::endl;
                printUsage(argv[0]);
                return 1;
            }

            std::string filename = argv[2];
            int source = std::stoi(argv[3]);
            int sink = std::stoi(argv[4]);
            int runs = (argc > 5) ? std::stoi(argv[5]) : 5;

            FlowNetwork network = FlowNetwork::loadFromFile(filename);
            std::vector<FlowNetwork> networks = {network};
            
            std::cout << "Note: Benchmarking not yet available" << std::endl;
        }
        else {
            std::cout << "Error: Unknown command " << command << std::endl;
            printUsage(argv[0]);
            return 1;
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}