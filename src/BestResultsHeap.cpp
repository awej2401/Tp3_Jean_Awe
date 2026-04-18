#include "BestResultsHeap.h"
#include <stdexcept>
#include <algorithm>

BestResultsHeap::BestResultsHeap(int capacity)
    : data_(nullptr), size_(0), capacity_(capacity) {
    if (capacity_ <= 0) {
        throw std::runtime_error("BestResultsHeap capacity must be > 0");
    }
    data_ = new AlignmentResult[capacity_];
}

BestResultsHeap::BestResultsHeap(BestResultsHeap&& other) noexcept
    : data_(other.data_), size_(other.size_), capacity_(other.capacity_) {
    other.data_ = nullptr;
    other.size_ = 0;
    other.capacity_ = 0;
}

BestResultsHeap& BestResultsHeap::operator=(BestResultsHeap&& other) noexcept {
    if (this != &other) {
        delete[] data_;

        data_ = other.data_;
        size_ = other.size_;
        capacity_ = other.capacity_;

        other.data_ = nullptr;
        other.size_ = 0;
        other.capacity_ = 0;
    }
    return *this;
}

BestResultsHeap::~BestResultsHeap() {
    delete[] data_;
}

bool BestResultsHeap::worseThan(const AlignmentResult& a,
                                const AlignmentResult& b) const {
    return a.mismatches > b.mismatches;
}

bool BestResultsHeap::betterThan(const AlignmentResult& a,
                                 const AlignmentResult& b) const {
    return a.mismatches < b.mismatches;
}

void BestResultsHeap::siftUp(int idx) {
    while (idx > 0) {
        int parent = (idx - 1) / 2;
        if (!worseThan(data_[idx], data_[parent])) {
            break;
        }
        std::swap(data_[idx], data_[parent]);
        idx = parent;
    }
}

void BestResultsHeap::siftDown(int idx) {
    while (true) {
        int left  = 2 * idx + 1;
        int right = 2 * idx + 2;
        int worst = idx;

        if (left < size_ && worseThan(data_[left], data_[worst])) {
            worst = left;
        }
        if (right < size_ && worseThan(data_[right], data_[worst])) {
            worst = right;
        }
        if (worst == idx) {
            break;
        }
        std::swap(data_[idx], data_[worst]);
        idx = worst;
    }
}

void BestResultsHeap::insert(const AlignmentResult& result) {
    if (size_ < capacity_) {
        data_[size_] = result;
        siftUp(size_);
        ++size_;
        return;
    }

    if (betterThan(result, data_[0])) {
        data_[0] = result;
        siftDown(0);
    }
}

AlignmentResult BestResultsHeap::getBest() const {
    if (size_ == 0) {
        throw std::runtime_error("BestResultsHeap::getBest on empty heap");
    }
    int bestIdx = 0;
    for (int i = 1; i < size_; ++i) {
        if (data_[i].mismatches < data_[bestIdx].mismatches) {
            bestIdx = i;
        }
    }
    return data_[bestIdx];
}

bool BestResultsHeap::isEmpty() const {
    return size_ == 0;
}

int BestResultsHeap::size() const {
    return size_;
}

int BestResultsHeap::capacity() const {
    return capacity_;
}
