#pragma once
#include "AlignmentResult.h"
#include <cstdint>

// Fixed-capacity max-heap ordered by mismatch count.
//
// We keep the N BEST (lowest-mismatch) results seen so far, where N is the
// capacity passed to the constructor (default 15).
// A max-heap: the root is always the WORST result
// AMONG those currently stored, so we can check in O(1)
// whether a new candidate is worth
// keeping, and evict the root in O(log N) when it is.
//
// 
// The underlying storage is a raw heap-allocated array of AlignmentResult.
//
class BestResultsHeap {
public:
    explicit BestResultsHeap(int capacity = 15);
    ~BestResultsHeap();

    void            insert(const AlignmentResult& result);
    AlignmentResult getBest() const;
    bool            isEmpty() const;
    int             size() const;
    int             capacity() const;

    // ❌ Interdire la copie
    BestResultsHeap(const BestResultsHeap&) = delete;
    BestResultsHeap& operator=(const BestResultsHeap&) = delete;

    // ✅ Autoriser le déplacement
    BestResultsHeap(BestResultsHeap&& other) noexcept;
    BestResultsHeap& operator=(BestResultsHeap&& other) noexcept;

private:
    AlignmentResult* data_;
    int              size_;
    int              capacity_;

    void siftUp(int idx);
    void siftDown(int idx);
    bool worseThan(const AlignmentResult& a, const AlignmentResult& b) const;
    bool betterThan(const AlignmentResult& a, const AlignmentResult& b) const;
};
