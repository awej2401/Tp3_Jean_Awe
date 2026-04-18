#include "Node.h"

Node::Node(uint64_t kmerCode_, const GenomicPosition& position_)
    : kmerCode(kmerCode_), position(position_), next(0) {}
