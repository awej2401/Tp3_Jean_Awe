#pragma once
#include "LinkedList.h"
#include "GenomicPosition.h"
#include "BaseLUT.h"
#include <string>
#include <cstdint>

// Hash table: kmer (2-bit encoded) -> linked list of (kmerCode, GenomicPosition).
// Bucket index  = kmerCode % capacity.
// Collision chains store the kmerCode in every Node so lookup() callers can
// skip nodes that belong to a different kmer sharing the same bucket.
//
// buckets_ is a raw heap-allocated array of LinkedList* (not a std::vector)
// so manage the allocation explicitly.
class HashTable {
public:
    explicit HashTable(uint32_t capacity);
    ~HashTable();

    uint8_t    encodeBase(char c) const;
    bool       encodeKmer(const std::string& kmer, uint64_t& code) const;
    void       insert(uint64_t kmerCode, const GenomicPosition& position);
    LinkedList* lookup(uint64_t kmerCode) const;
    uint32_t   capacity() const;

    HashTable(const HashTable&) = delete;
    HashTable& operator=(const HashTable&) = delete;

    HashTable(HashTable&& other) noexcept;
    HashTable& operator=(HashTable&& other) noexcept;

private:
    LinkedList** buckets_;
    uint32_t     capacity_;
};