#include "sequential.h"
#include <limits>
#include <algorithm>
#include <iomanip>
#include <vector>



int PushRelabelSequential::maxFlow(FlowNetwork& network, int source, int sink) {
    int n = network.getNumVertices();
    std::cout << "Starting Sequential Push-Relabel with " << n
              << " vertices (source=" << source << ", sink=" << sink << ")" << std::endl;

    if (source < 0 || source >= n || sink < 0 || sink >= n || source == sink) {
        throw std::invalid_argument("Invalid source or sink");
    }

    auto& graph = const_cast<std::vector<std::vector<FlowNetwork::Edge>>&>(network.getGraph());
    std::vector<int> excess(n, 0), height(n, 0);
    std::queue<int> active_vertices;
    std::vector<bool> in_queue(n, false);

    int relabel_since_last = 0;
    const int global_freq = n;
    int total_global = 0;

    initialize(graph, excess, height, active_vertices, in_queue, source, sink, n);

    globalRelabel(graph, height, source, sink, n);
    total_global++;
    active_vertices = {};
    std::fill(in_queue.begin(), in_queue.end(), false);
    for (int i = 0; i < n; ++i) {
        if (i != source && i != sink && excess[i] > 0 && height[i] < n) {
            active_vertices.push(i);
            in_queue[i] = true;
        }
    }

    int iterations = 0, total_relabels = 0;
    while (!active_vertices.empty() && iterations < MAX_ITERATIONS) {
        if (relabel_since_last >= global_freq) {
            globalRelabel(graph, height, source, sink, n);
            total_global++;
            relabel_since_last = 0;

            std::queue<int> new_queue;
            std::fill(in_queue.begin(), in_queue.end(), false);
            while (!active_vertices.empty()) {
                int u = active_vertices.front(); active_vertices.pop();
                if (u != source && u != sink && excess[u] > 0 && height[u] < n) {
                    new_queue.push(u);
                    in_queue[u] = true;
                }
            }
            for (int i = 0; i < n; ++i) {
                if (i != source && i != sink && excess[i] > 0 && height[i] < n && !in_queue[i]) {
                    new_queue.push(i);
                    in_queue[i] = true;
                }
            }
            active_vertices = std::move(new_queue);
            if (active_vertices.empty()) break;
        }

        int u = active_vertices.front();
        active_vertices.pop();
        in_queue[u] = false;

        int relabels = discharge(graph, excess, height, active_vertices, in_queue, u, source, sink, n);
        total_relabels += relabels;
        relabel_since_last += relabels;
        iterations++;
    }

    if (iterations >= MAX_ITERATIONS) {
        std::cout << "WARNING: Reached maximum iterations. Max flow might be incorrect." << std::endl;
    }

    int max_flow = excess[sink];
    std::cout << "Completed in " << iterations << " iterations. Total local relabels: "
              << total_relabels << ", global relabels: " << total_global
              << ". Max flow: " << max_flow << std::endl;
    return max_flow;
}

void PushRelabelSequential::initialize(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                     std::vector<int>& excess, std::vector<int>& height,
                                     std::queue<int>& active_vertices, std::vector<bool>& in_queue,
                                     int source, int sink, int n) {
    height.assign(n, 0);
    excess.assign(n, 0);
    in_queue.assign(n, false);
    height[source] = n;

    for (auto& edge : const_cast<std::vector<FlowNetwork::Edge>&>(graph[source])) {
        if (edge.capacity > 0) {
            int flow = edge.capacity;
            int v = edge.to;
            auto& rev = const_cast<std::vector<FlowNetwork::Edge>&>(graph[v])[edge.rev];
            edge.flow += flow;
            rev.flow -= flow;
            excess[v] += flow;
            excess[source] -= flow;
            if (v != source && v != sink && !in_queue[v]) {
                active_vertices.push(v);
                in_queue[v] = true;
            }
        }
    }
}

bool PushRelabelSequential::push(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                               std::vector<int>& excess, std::vector<int>& height,
                               std::queue<int>& active_vertices, std::vector<bool>& in_queue,
                               int u, int idx, int source, int sink) {
    auto& e = graph[u][idx];
    int v = e.to;
    if (e.capacity - e.flow <= 0 || height[u] != height[v] + 1) return false;
    auto& rev = graph[v][e.rev];
    int delta = std::min(excess[u], e.capacity - e.flow);
    if (delta <= 0) return false;
    e.flow += delta;
    rev.flow -= delta;
    excess[u] -= delta;
    excess[v] += delta;
    if (v != source && v != sink && excess[v] > 0 && !in_queue[v]) {
        active_vertices.push(v);
        in_queue[v] = true;
    }
    return true;
}

bool PushRelabelSequential::relabel(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                  std::vector<int>& height, int u, int n) {
    int min_h = INF_HEIGHT;
    for (auto& e : graph[u]) {
        if (e.capacity - e.flow > 0)
            min_h = std::min(min_h, height[e.to]);
    }
    height[u] = (min_h < INF_HEIGHT ? min_h + 1 : n);
    return true;
}

int PushRelabelSequential::discharge(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                   std::vector<int>& excess, std::vector<int>& height,
                                   std::queue<int>& active_vertices, std::vector<bool>& in_queue,
                                   int u, int source, int sink, int n) {
    if (u == source || u == sink || excess[u] <= 0) return 0;
    int relabels = 0;
    while (excess[u] > 0 && height[u] < n) {
        bool did_push = false;
        for (size_t i = 0; i < graph[u].size() && excess[u] > 0; ++i) {
            if (push(graph, excess, height, active_vertices, in_queue, u, i, source, sink))
                did_push = true;
        }
        if (!did_push && excess[u] > 0) {
            relabel(graph, height, u, n);
            relabels++;
        }
    }
    return relabels;
}

void PushRelabelSequential::globalRelabel(std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                         std::vector<int>& height, int source, int sink, int n) {
    height.assign(n, n);
    height[sink] = 0;
    std::queue<int> q;
    q.push(sink);
    std::vector<bool> vis(n, false);
    vis[sink] = true;
    while (!q.empty()) {
        int v = q.front(); q.pop();
        for (auto& e : graph[v]) {
            int u = e.to;
            int rev = e.rev;
            if (u >= 0 && u < n && rev < graph[u].size()) {
                auto& rev_e = graph[u][rev];
                if (rev_e.capacity - rev_e.flow > 0 && !vis[u]) {
                    vis[u] = true;
                    height[u] = height[v] + 1;
                    q.push(u);
                }
            }
        }
    }
    height[source] = n;
}

void PushRelabelSequential::printState(const std::vector<std::vector<FlowNetwork::Edge>>& graph,
                                     const std::vector<int>& excess, const std::vector<int>& height,
                                     const std::queue<int>& active_vertices, int n) {
    std::cout << "--- Current State ---" << std::endl;
    std::cout << "  Heights: ";
    for (int i = 0; i < n; i++) {
        std::cout << (height[i] >= n ? "INF " : std::to_string(height[i]) + " ");
    }
    std::cout << std::endl;
    if (!excess.empty()) {
        std::cout << "  Excess:  ";
        for (int i = 0; i < n; i++) std::cout << excess[i] << " ";
        std::cout << std::endl;
    }
    std::cout << "  Active Queue (" << active_vertices.size() << "): ";
    std::queue<int> tmp = active_vertices;
    while (!tmp.empty()) { std::cout << tmp.front() << " "; tmp.pop(); }
    std::cout << std::endl;
    std::cout << "---------------------" << std::endl;
}
