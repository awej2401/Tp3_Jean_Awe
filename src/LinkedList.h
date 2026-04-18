#pragma once
#include "Node.h"
#include "GenomicPosition.h"
#include <cstdint>

// Singly-linked list of Nodes used as a hash-bucket collision chain.
class LinkedList {
public:
    LinkedList();
    ~LinkedList();

    void insert(uint64_t kmerCode, const GenomicPosition& position);
    Node* first() const;

    LinkedList(const LinkedList&) = delete;
    LinkedList& operator=(const LinkedList&) = delete;

    LinkedList(LinkedList&& other) noexcept;
    LinkedList& operator=(LinkedList&& other) noexcept;

private:
    Node* head_;
};
