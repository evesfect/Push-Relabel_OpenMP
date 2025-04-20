# Push-Relabel_OpenMP


Design choices for the sequential Push-Relabel

FIFO vs Highest Label

Time Complexity
    FIFO: O(V^3), HL: O((V^2)E^(1/2))

Vertex Selection Data Structure
    FIFO: Queue O(1) access, HL: Priority Queue O(log n) access / Bucket

Design choices for the sequential Push-Relabel

Global Relabeling with Backwards Parallel Breadth-First Search

Batching with thread concurrent queue

