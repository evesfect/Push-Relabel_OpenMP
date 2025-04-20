#include "test.h"
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>

// Convert flow network to Boost graph
boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
                     boost::no_property,
                     boost::property<boost::edge_capacity_t, int,
                     boost::property<boost::edge_residual_capacity_t, int,
                     boost::property<boost::edge_reverse_t, boost::graph_traits<
                         boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
                                             boost::no_property,
                                             boost::property<boost::edge_capacity_t, int,
                                             boost::property<boost::edge_residual_capacity_t, int,
                                             boost::property<boost::edge_reverse_t, void*>
                                             >>
                                            >
                         >::edge_descriptor
                     >>>>

convertToBoostGraph(const FlowNetwork& network) {
    using namespace boost;
    
    typedef adjacency_list<vecS, vecS, directedS,
                         no_property,
                         property<edge_capacity_t, int,
                         property<edge_residual_capacity_t, int,
                         property<edge_reverse_t, 
                                 graph_traits<adjacency_list<vecS, vecS, directedS,
                                                           no_property,
                                                           property<edge_capacity_t, int,
                                                           property<edge_residual_capacity_t, int,
                                                           property<edge_reverse_t, void*>
                                                           >>
                                                          >
                                 >::edge_descriptor
                         >>>> Graph;
    
    Graph g(network.getNumVertices());
    
    property_map<Graph, edge_capacity_t>::type capacity = get(edge_capacity, g);
    property_map<Graph, edge_residual_capacity_t>::type residual_capacity = get(edge_residual_capacity, g);
    property_map<Graph, edge_reverse_t>::type rev = get(edge_reverse, g);
    
    const auto& graph = network.getGraph();
    
    // Add edges to the boost graph
    for (int from = 0; from < network.getNumVertices(); ++from) {
        for (const auto& edge : graph[from]) {
            if (edge.capacity > 0) { // Only add forward edges
                graph_traits<Graph>::edge_descriptor e1, e2;
                bool inserted;
                
                // Add the forward edge
                tie(e1, inserted) = add_edge(from, edge.to, g);
                capacity[e1] = edge.capacity;
                
                // Add the reverse edge for residual capacity
                tie(e2, inserted) = add_edge(edge.to, from, g);
                capacity[e2] = 0;
                
                // Store reverse edge pointers
                rev[e1] = e2;
                rev[e2] = e1;
            }
        }
    }
    
    return g;
}

int MaxFlowTester::computeExpectedFlow(const FlowNetwork& network, int source, int sink) {
    using namespace boost;
    auto g = convertToBoostGraph(network);
    
    return push_relabel_max_flow(g, source, sink);
}

bool MaxFlowTester::verifyCorrectness(const FlowNetwork& network, int source, int sink, int computed_flow) {
    int expected_flow = computeExpectedFlow(network, source, sink);
    
    if (computed_flow == expected_flow) {
        std::cout << "Test passed: Max flow = " << computed_flow << std::endl;
        return true;
    } else {
        std::cout << "Test failed: Computed = " << computed_flow 
                  << ", Expected = " << expected_flow << std::endl;
        return false;
    }
}

// TODO:
bool MaxFlowTester::compareImplementations(const FlowNetwork& network, int source, int sink) {
    std::cout << "Implementation comparison not yet available." << std::endl;
    return false;
}

// TODO:
void MaxFlowTester::runBenchmark(const std::vector<FlowNetwork>& networks, 
                               int source, int sink,
                               int num_runs) {
    std::cout << "Benchmarking not yet available." << std::endl;
}