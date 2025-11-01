# Optimizing and Parallelizing the Push-Relabel Algorithm

**Ege Yurtsever** | University of Bergen | 2025

A custom implementation and parallelization of the Goldberg-Tarjan Push-Relabel maximum flow algorithm using OpenMP. This project demonstrates that Global Relabeling heuristics provide ~46x sequential speedup over basic FIFO, and parallelization provides an additional 2-4x speedup on large dense graphs. Inspired by prior work on parallel push-relabel algorithms (e.g., Baumstark et al.), this implementation explores an alternative approach that showed superior performance in my experiments with Intel Xeon and Core i7 processors.

## What's Implemented

Three versions of Push-Relabel with different optimization levels:

1.  **Sequential FIFO with Global Relabeling** (`src/sequential.cpp`)
    
    -   `PushRelabelSequential::maxFlow()`
    -   Backward BFS-based Global Relabeling every n iterations
    -   ~46x faster than basic FIFO
2.  **Parallel Global Relabeling, Sequential Discharge** (`src/parallel.cpp`)
    
    -   `PushRelabelParallel::maxFlow()`
    -   OpenMP-parallelized level-synchronous BFS
    -   Hybrid approach: sequential for n<1300, parallel otherwise
3.  **Fully Parallel** (`src/parallel.cpp`)
    
    -   `PushRelabelParallel::maxFlow_ActiveParallel()`
    -   Atomic operations on excess flow
    -   Batch-parallel discharge with thread-local active lists
    -   Best for large dense graphs (>20k vertices, density >0.3)

## Performance Results

### Benchmark Data (8 threads, edge density = 0.3)

(FIFO becomes untestable for the graph sizes above 5000 due to incredible time cost)(GR = Global Relabeling Heuristics, DC = Discharge)

| Graph Size | FIFO (ms) | FIFO+GR (ms) | Parallel GR (ms) | Parallel DC (ms) |
|------------|------------|----------------|------------------|-----------------|
| 1000       | 6          | 5              | 5                | 7               |
| 2000       | 2290       | 43             | 32               | 42              |
| 5000       | 38109      | 305            | 168              | 164             |
| 10000      | -          | 793            | 295              | 320             |
| 20000      | -          | 8289           | 4842             | 3941            |
| 30000      | -          | 13885          | 8762             | 6851            |
| 40000      | -          | 63677          | 23164            | 16459           |



![alt text](image.png)

### Key Observations

-   **Global Relabeling is critical**: 46x speedup over basic FIFO (tested on graphs with n=1000-2000, d=0.1-0.5)
-   **Cost breakdown (n=10000, d=0.3)**: Global Relabels 79.73%, Push 3.37%, Relabel 2.94%, Queue 0.62%
-   **Parallelizing Global Relabeling**: 2-3x speedup on large graphs
-   **Full parallelization**: Additional speedup only on dense graphs; overhead dominates on sparse graphs
-   **Sparse graphs (d=0.01)**: Full parallelization unreliable, use Parallel GR only
-   **Dense graphs (d=0.3)**: Full parallelization scales well with size

## Setup

**Requirements:**

-   C++17 compiler with OpenMP support
-   Boost Graph Library (for correctness testing against Boost's push-relabel)

**Configure Boost:**Set `BOOST_ROOT` environment variable or edit `Makefile` line 13-17 if Boost is not in default location.

## Build

```bash
make parallel              # Main executable: bin/max_flowmake benchmark-parallel    # Performance tests: bin/performance_testmake test                  # Correctness tests (runs automatically)make clean
```

**Note:** Only these targets are supported. Do not use `make sequential`.

## Usage

### Generate and Test Graphs

```bash
# Generate random graph
bin/max_flow -g 1000 0.3 1 100 bin/graph.dat
# Verify with Boost
bin/max_flow -t bin/graph.dat 0 999
```

### Run Benchmarks

```bash
# Default: 1000 vertices, density 0.3, 8 threads
bin/performance_test
# Custom
bin/performance_test 20000 0.3 8
```

Benchmark output compares all three implementations (Sequential+GR, Parallel GR, Fully Parallel).

## Implementation Details

### Core Data Structures

**FlowNetwork** (`src/flow_network.h`)

-   Adjacency list with forward/reverse edge pairs
-   Each edge stores: `to`, `capacity`, `flow`, `rev` (reverse edge index)
-   `graph[u][i].rev` points to reverse edge index in `graph[v]`

**Algorithm State** (Sequential and Parallel)

-   `excess[]`: Excess flow at each vertex
-   `height[]`: Distance labels (updated by Global Relabeling)
-   `active_queue`: FIFO queue of vertices with positive excess
-   `in_queue[]`: Tracks queue membership

### Global Relabeling

**Sequential** (`PushRelabelSequential::globalRelabel` in `src/sequential.cpp`):

-   Backward BFS from sink in residual graph
-   Sets height[v] = shortest path distance to sink
-   Triggered every n iterations (n = |V|)

**Parallel** (`PushRelabelParallel::globalRelabel` in `src/parallel.cpp`):

-   Level-synchronous BFS with OpenMP
-   Each thread processes portion of current frontier
-   Thread-local next frontiers merged after each level
-   Critical section for height updates (uses `height[u] == n` check to avoid duplicates)
-   Hybrid: sequential for n<1300, parallel otherwise

### Discharge Operations

**Sequential Discharge** (`PushRelabelSequential::discharge` in `src/sequential.cpp`):

-   While vertex has excess flow, try to push to neighbors at height[u] - 1
-   If no valid push exists, relabel (increase height[u])

**Parallel Discharge** (`PushRelabelParallel::discharge_Atomic` in `src/parallel.cpp`):

-   Batch-parallel approach: extract batch from active_queue
-   Each thread discharges subset with atomic excess updates
-   Thread-local newly-active lists merged back to global queue
-   Re-check batch vertices for remaining excess

Key difference: Atomic operations on `excess[]` allow concurrent pushes without data races. Critical sections only for queue management.

### Key Challenges in Parallelization

1.  **Atomic contention**: `std::atomic<int>` on excess values causes performance degradation on sparse graphs
2.  **Queue synchronization**: Active queue updates require critical sections
3.  **Work distribution**: Unbalanced active vertex distribution affects parallel efficiency
4.  **Memory access patterns**: Random access in discharge operations limits cache effectiveness

These issues explain why full parallelization only helps on large dense graphs where computation dominates overhead.

## References

**Technical Report:** `Parallelizing and Optimizing the Push-Relabel Algorithm-University of Bergen-Ege Yurtsever.pdf`

**Papers:**

-   Cormen et al., *Introduction to Algorithms*, 2nd ed., MIT Press, 2001
-   Baumstark, Blelloch, Shun, *Efficient Implementation of a Synchronous Parallel Push-Relabel Algorithm*, CMU/KIT, 2001
