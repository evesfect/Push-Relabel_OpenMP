#include "flow_network.h"

FlowNetwork::FlowNetwork(int n) : num_vertices(n), graph(n) {}

void FlowNetwork::addEdge(int from, int to, int capacity) {
    // Add forward edge
    graph[from].push_back(Edge(to, capacity, 0, graph[to].size()));
    
    // Add residual edge w/ 0 capacity
    graph[to].push_back(Edge(from, 0, 0, graph[from].size() - 1));
}

void FlowNetwork::resetFlow() {
    for (auto& edges : graph) {
        for (auto& edge : edges) {
            edge.flow = 0;
        }
    }
}

int FlowNetwork::getNumVertices() const {
    return num_vertices;
}

const std::vector<std::vector<FlowNetwork::Edge>>& FlowNetwork::getGraph() const {
    return graph;
}

void FlowNetwork::print() const {
    std::cout << "Flow Network with " << num_vertices << " vertices:" << std::endl;
    for (int from = 0; from < num_vertices; ++from) {
        for (const auto& edge : graph[from]) {
            if (edge.capacity > 0) { // Only print forward edges, not residual
                std::cout << from << " -> " << edge.to << ": capacity = " 
                          << edge.capacity << ", flow = " << edge.flow << std::endl;
            }
        }
    }
}

bool FlowNetwork::saveToFile(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file) {
        return false;
    }

    file << num_vertices << std::endl;
    
    // Write forward edges
    for (int from = 0; from < num_vertices; ++from) {
        for (const auto& edge : graph[from]) {
            if (edge.capacity > 0) {
                file << from << " " << edge.to << " " << edge.capacity << std::endl;
            }
        }
    }
    
    file.close();
    return true;
}

FlowNetwork FlowNetwork::loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    int n;
    file >> n;
    FlowNetwork network(n);

    int from, to, capacity;
    while (file >> from >> to >> capacity) {
        network.addEdge(from, to, capacity);
    }

    file.close();
    return network;
}