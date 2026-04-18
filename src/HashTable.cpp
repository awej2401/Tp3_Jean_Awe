#include "HashTable.h"
#include <stdexcept>

HashTable::HashTable(uint32_t capacity)
    : buckets_(nullptr), capacity_(capacity) {
    if (capacity_ == 0u) {
        throw std::runtime_error("HashTable capacity must be > 0");
    }
    buckets_ = new LinkedList*[capacity_];
    for (uint32_t i = 0; i < capacity_; ++i) {
        buckets_[i] = new LinkedList();
    }
}

HashTable::HashTable(HashTable&& other) noexcept
    : buckets_(other.buckets_), capacity_(other.capacity_) {
    other.buckets_ = nullptr;
    other.capacity_ = 0;
}

HashTable& HashTable::operator=(HashTable&& other) noexcept {
    if (this != &other) {
        if (buckets_ != nullptr) {
            for (uint32_t i = 0; i < capacity_; ++i) {
                delete buckets_[i];
            }
            delete[] buckets_;
        }

        buckets_ = other.buckets_;
        capacity_ = other.capacity_;

        other.buckets_ = nullptr;
        other.capacity_ = 0;
    }
    return *this;
}

HashTable::~HashTable() {
    if (buckets_ != nullptr) {
        for (uint32_t i = 0; i < capacity_; ++i) {
            delete buckets_[i];
        }
        delete [] buckets_;
    }
}

uint8_t HashTable::encodeBase(char c) const {
    return BASE_LUT[static_cast<unsigned char>(c)];
}

bool HashTable::encodeKmer(const std::string& kmer, uint64_t& code) const {
    code = 0u;
    for (uint32_t i = 0; i < static_cast<uint32_t>(kmer.size()); ++i) {
        uint8_t b = encodeBase(kmer[i]);
        if (b == INVALID) {
            code = 0u;
            return false;
        }
        code = (code << 2) | static_cast<uint64_t>(b);
    }
    return true;
}

void HashTable::insert(uint64_t kmerCode, const GenomicPosition& position) {
    uint32_t idx = static_cast<uint32_t>(kmerCode % capacity_);
    buckets_[idx]->insert(kmerCode, position);
}

LinkedList* HashTable::lookup(uint64_t kmerCode) const {
    uint32_t idx = static_cast<uint32_t>(kmerCode % capacity_);
    return buckets_[idx];
}

uint32_t HashTable::capacity() const {
    return capacity_;
}
