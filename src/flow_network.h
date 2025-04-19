#ifndef FLOW_NETWORK_H
#define FLOW_NETWORK_H

#include <vector>
#include <iostream>
#include <fstream>
#include <string>

class FlowNetwork {
public:
    struct Edge {
        int to;
        int capacity;
        int flow;
        int rev; // Index of reverse edge

        Edge(int t, int c, int f, int r) : to(t), capacity(c), flow(f), rev(r) {}
    };

    FlowNetwork(int n);
    void addEdge(int from, int to, int capacity);
    void resetFlow();
    int getNumVertices() const;
    const std::vector<std::vector<Edge>>& getGraph() const;

    // For testing and visualization
    void print() const;
    bool saveToFile(const std::string& filename) const;
    static FlowNetwork loadFromFile(const std::string& filename);

private:
    int num_vertices;
    std::vector<std::vector<Edge>> graph;
};

#endif // FLOW_NETWORK_H