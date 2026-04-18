#pragma once
#include "GenomicPosition.h"
#include <cstdint>

// Singly-linked list node.
// Stores the 2-bit encoded kmer code alongside the GenomicPosition so that
// hash-bucket chain walkers can reject nodes from colliding kmers.
class Node {
public:
    uint64_t        kmerCode;
    GenomicPosition position;
    Node*           next;

    Node(uint64_t kmerCode, const GenomicPosition& position);
};
