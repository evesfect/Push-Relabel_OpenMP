# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an implementation of the Push-Relabel maximum flow algorithm with parallel optimizations using OpenMP. The project compares sequential and parallel implementations of the algorithm, with a focus on performance improvements through global relabeling heuristics and parallelization strategies.

**Key Research Note**: Global Relabeling heuristic provides ~45x speedup for sequential and ~200x speedup for parallel compared to basic FIFO. This is the critical optimization that makes the algorithm practical.

## Build Commands

**Prerequisites**: Boost library must be installed. Set `BOOST_ROOT` environment variable or adjust the Boost directory path in the Makefile.

### Primary Build Targets

```bash
# Build parallel version (primary use case)
make parallel

# Build for performance benchmarking (recommended for testing)
make benchmark-parallel

# Run correctness tests
make test

# Clean build artifacts
make clean
```

**Important**: `make sequential` and similar commands are not supported. Only use the commands listed above.

### Running the Executables

After building with `make parallel`, the main executable is located at `bin/max_flow`:

```bash
# Generate a random graph
bin/max_flow -g <vertices> <edge_probability> <min_capacity> <max_capacity> <output_file>
# Example: bin/max_flow -g 1000 0.3 1 100 bin/large_graph.dat

# Test max flow on a graph file (computes expected flow using Boost)
bin/max_flow -t <graph_file> <source> <sink>

# Benchmark (not fully implemented)
bin/max_flow -b <graph_file> <source> <sink> [runs]
```

After building with `make benchmark-parallel`:

```bash
# Run performance tests with custom parameters
bin/performance_test [vertices] [density] [num_threads]
# Example: bin/performance_test 5000 0.3 8
# Defaults: 1000 vertices, 0.3 density, 8 threads
```

After building with `make test`:

```bash
# Runs correctness tests automatically
```

## Architecture

### Core Components

**FlowNetwork** (`src/flow_network.h/cpp`):
- Residual graph representation using adjacency lists
- Each edge stores: destination vertex, capacity, current flow, and reverse edge index
- Methods for graph construction, flow reset, and file I/O

**PushRelabelSequential** (`src/sequential.h/cpp`):
- FIFO-based sequential implementation
- Uses Global Relabeling heuristic (critical for performance)
- Key operations: `initialize()`, `push()`, `relabel()`, `discharge()`

**PushRelabelParallel** (`src/parallel.h/cpp`):
- Two parallel strategies:
  1. `maxFlow()`: Original parallel approach with standard parallelization
  2. `maxFlow_ActiveParallel()`: Active vertex parallelization with atomic operations on excess values
- Parallel Global Relabeling using BFS with OpenMP
- Thread count configurable (defaults to system max if not specified)

**FlowNetworkGenerator** (`src/generator.h/cpp`):
- Random graph generation with configurable parameters
- Used for testing and benchmarking

**MaxFlowTester** (`src/test.h/cpp`):
- Verification against Boost's push-relabel implementation
- Validates flow conservation and capacity constraints

### Algorithm Flow

1. **Initialization**: Saturate all edges from source, set height[source] = n
2. **Main Loop**: Process active vertices (vertices with excess flow) using discharge operations
3. **Discharge**: Attempt to push excess to lower neighbors; relabel if no valid push exists
4. **Global Relabeling**: Periodically perform BFS from sink to update height labels (crucial optimization)
5. **Termination**: When no active vertices remain, flow at sink is maximum flow

### Parallelization Strategy

The `maxFlow_ActiveParallel()` method uses:
- Atomic operations on `excess[]` array to handle concurrent updates
- Thread-local lists for newly activated vertices (merged after parallel section)
- Critical sections for updating the global active queue
- Parallel global relabeling using OpenMP-parallelized BFS

**Performance Note**: According to development notes, simple parallelization shows good speedup only on large graphs due to atomic/contention overhead. Most time is spent in push operations.

### Testing Structure

**Correctness Tests** (`test/correctness_test.cpp`):
- `testTinyGraph()`: 4-vertex manual test case
- `testSmallGraph()`: 6-vertex manual test case
- `testRandomGraph()`: 1000-vertex random graph
- Validates against Boost's implementation using `MaxFlowTester`

**Performance Tests** (`test/performance_test.cpp`):
- `testGlobalRelabelingPerformance()`: Isolates global relabeling performance
- `testMaxFlowPerformance()`: Compares Sequential, Parallel (Original), and Parallel (ActiveParallel)
- Uses `Timer` class for measurements (defined in `src/test.h`)
- Configurable via command-line arguments

## Development Guidelines

### Compilation Flags

- `-std=c++17`: C++17 standard
- `-O3`: Aggressive optimization
- `-fopenmp`: OpenMP support for parallel versions
- Boost include path configured via `BOOST_INC` variable

### When Modifying Algorithms

- Maintain consistency between sequential and parallel versions for correctness comparison
- Global Relabeling is essential for performance—do not remove it
- Be cautious with atomic operations; excessive use causes performance degradation
- Test with both small manual graphs and large random graphs
- Verify correctness against Boost implementation before performance tuning

### File Naming Convention

The repository uses inconsistent casing: `MakeFile` exists but should be `Makefile`. When referencing the build file, use `make` commands which work regardless of the filename casing on case-insensitive systems.
