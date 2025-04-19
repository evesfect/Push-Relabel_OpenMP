#ifndef SEQUENTIAL_H
#define SEQUENTIAL_H

#include "flow_network.h"

class PushRelabelSequential {
public:
    // Find maximum flow using push-relabel algorithm
    static int maxFlow(FlowNetwork& network, int source, int sink);
};

#endif // SEQUENTIAL_H