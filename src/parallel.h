#ifndef PARALLEL_H
#define PARALLEL_H

#include "flow_network.h"

// Only include OpenMP header if OpenMP is enabled
#ifdef _OPENMP
#include <omp.h>
#endif

class PushRelabelParallel {
public:
    // Find maximum flow using parallel push-relabel algorithm
    static int maxFlow(FlowNetwork& network, int source, int sink, int num_threads = 0);
};

#endif // PARALLEL_H