/*
 * main.cpp
 *
 * Compile:  make
 * Usage:    ./bin/genome_mapper genome.fasta reads.fasta [OPTIONS]
 *
 * Required:
 *   genome.fasta          genome reference file
 *   reads.fasta           reads to align
 *
 * Options:
 *   --kmer_size   INT     k-mer length used for indexing (default: 21, max: 31)
 *   --table_size  INT     number of hash-table buckets  (default: 4000000)
 *   --heap_size   INT     max candidates per read kept in the BestResultsHeap (default: 15)
 *   --help                print this help message and exit
 *
 * Output columns (TSV, stdout):
 *   #read  chrom  start(0-based)  end(0-based inclusive)  mismatches  matches  candidates
 *   Unmapped reads: chrom=*  start=end=mismatches=matches=-1  candidates=0
 */

#include "BaseLUT.h"
#include "HashTable.h"
#include "GenomicPosition.h"
#include "AlignmentResult.h"
#include "BestResultsHeap.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <cstdlib>
#include <cstdint>
#include <set>
#include <utility>

static void sanitizeSeq(std::string& seq) {
    for (uint32_t i = 0; i < static_cast<uint32_t>(seq.size()); ++i) {
        unsigned char c = static_cast<unsigned char>(seq[i]);
        if (BASE_LUT[c] == INVALID)
            seq[i] = 'X';
    }
}

typedef std::vector<std::pair<std::string, std::string> > GenomeMap;

static GenomeMap readFasta(const std::string& path) {
    std::ifstream in(path.c_str());
    if (!in.is_open())
        throw std::runtime_error("Cannot open file: " + path);

    GenomeMap records;
    std::string line, name, seq;

    while (std::getline(in, line)) {
        if (!line.empty() && line[line.size() - 1] == '\r')
            line.erase(line.size() - 1);
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!name.empty()) {
                sanitizeSeq(seq);
                records.push_back(std::make_pair(name, seq));
            }
            name = line.substr(1);
            uint32_t sp = static_cast<uint32_t>(name.find_first_of(" \t"));
            if (sp != static_cast<uint32_t>(std::string::npos))
                name = name.substr(0, sp);
            seq.clear();
        } else {
            std::transform(line.begin(), line.end(), line.begin(), ::toupper);
            seq += line;
        }
    }
    if (!name.empty()) {
        sanitizeSeq(seq);
        records.push_back(std::make_pair(name, seq));
    }
    return records;
}

static const std::string* findChrom(const GenomeMap& genome,
                                    const std::string& name) {
    for (uint32_t i = 0; i < static_cast<uint32_t>(genome.size()); ++i)
        if (genome[i].first == name) return &genome[i].second;
    return 0;
}

static void scoreAlignment(const std::string& read,
                           const std::string& chromSeq,
                           long long          genomeStart,
                           int&               mismatches,
                           int&               matches) {
    mismatches = 0;
    matches    = 0;
    const long long glen = static_cast<long long>(chromSeq.size());
    const int       rlen = static_cast<int>(read.size());

    for (int i = 0; i < rlen; ++i) {
        long long gpos = genomeStart + i;
        if (gpos < 0 || gpos >= glen) {
            ++mismatches;
            continue;
        }
        uint8_t r = BASE_LUT[static_cast<unsigned char>(read[i])];
        uint8_t g = BASE_LUT[static_cast<unsigned char>(chromSeq[static_cast<uint32_t>(gpos)])];
        if (r == INVALID || g == INVALID) continue;
        if (r == g) ++matches;
        else        ++mismatches;
    }
}

static void buildIndex(HashTable&       ht,
                       const GenomeMap& genome,
                       int              k) {
    uint32_t inserted = 0, skipped = 0;
    const uint64_t mask = (k < 32)
                          ? ((uint64_t(1) << (2 * k)) - 1)
                          : ~uint64_t(0);

    for (uint32_t ri = 0; ri < static_cast<uint32_t>(genome.size()); ++ri) {
        const std::string& chrom = genome[ri].first;
        const std::string& seq   = genome[ri].second;
        const uint32_t     len   = static_cast<uint32_t>(seq.size());
        if (len < static_cast<uint32_t>(k)) continue;

        uint64_t code = 0u;
        int      validRun = 0;
        for (uint32_t j = 0; j < len; ++j) {
            uint8_t b = BASE_LUT[static_cast<unsigned char>(seq[j])];
            if (b == INVALID) {
                code = 0u;
                validRun = 0;
                ++skipped;
                continue;
            }
            code = ((code << 2) | static_cast<uint64_t>(b)) & mask;
            ++validRun;
            if (validRun >= k) {
                uint32_t start = j + 1u - static_cast<uint32_t>(k);
                ht.insert(code, GenomicPosition(chrom, start));
                ++inserted;
            }
        }
    }
    std::cerr << "[index] k=" << k
              << "  inserted=" << inserted
              << "  skipped(non-ACGT)=" << skipped << "\n";
}

static void mapReads(const HashTable&  ht,
                     const GenomeMap&  genome,
                     const GenomeMap&  reads,
                     int               k,
                     int               heapSize) {
    uint32_t mapped = 0, unmapped = 0;
    std::cout << "#read\tchrom\tstart\tend\tmismatches\tmatches\tcandidates\n";

    for (uint32_t ri = 0; ri < static_cast<uint32_t>(reads.size()); ++ri) {
        const std::string& readName = reads[ri].first;
        const std::string& seq      = reads[ri].second;
        const int          rlen     = static_cast<int>(seq.size());

        BestResultsHeap heap(heapSize);
        std::set<std::pair<std::string, long long> > seen;

        if (rlen >= k) {
            const uint64_t mask = (k < 32)
                                  ? ((uint64_t(1) << (2 * k)) - 1)
                                  : ~uint64_t(0);
            uint64_t code = 0u;
            int validRun = 0;

            for (int i = 0; i < rlen; ++i) {
                uint8_t b = BASE_LUT[static_cast<unsigned char>(seq[i])];
                if (b == INVALID) {
                    code = 0u;
                    validRun = 0;
                    continue;
                }
                code = ((code << 2) | static_cast<uint64_t>(b)) & mask;
                ++validRun;
                if (validRun < k) continue;

                const int seedStart = i - k + 1;
                LinkedList* bucket = ht.lookup(code);
                for (Node* cur = bucket->first(); cur != 0; cur = cur->next) {
                    if (cur->kmerCode != code) continue;
                    const GenomicPosition& gp = cur->position;
                    long long readStart = static_cast<long long>(gp.offset) -
                                          static_cast<long long>(seedStart);
                    std::pair<std::string, long long> key(gp.chrom, readStart);
                    if (seen.find(key) != seen.end()) continue;
                    seen.insert(key);

                    if (readStart < 0) continue;
                    const std::string* chromSeq = findChrom(genome, gp.chrom);
                    if (chromSeq == 0) continue;

                    int mm = 0;
                    int m  = 0;
                    scoreAlignment(seq, *chromSeq, readStart, mm, m);
                    heap.insert(AlignmentResult(
                        GenomicPosition(gp.chrom, static_cast<uint32_t>(readStart)),
                        mm, m));
                }
            }
        }

        if (!heap.isEmpty()) {
            AlignmentResult best = heap.getBest();
            long long startPos = static_cast<long long>(best.position.offset);
            long long endPos   = startPos + static_cast<long long>(rlen) - 1;
            std::cout << readName
                      << "\t" << best.position.chrom
                      << "\t" << startPos
                      << "\t" << endPos
                      << "\t" << best.mismatches
                      << "\t" << best.matches
                      << "\t" << heap.size()
                      << "\n";
            ++mapped;
        } else {
            std::cout << readName << "\t*\t-1\t-1\t-1\t-1\t0\n";
            ++unmapped;
        }
    }

    std::cerr << "[map] mapped=" << mapped
              << "  unmapped=" << unmapped << "\n";
}

static void printHelp(const char* progName) {
    std::cout
        << "\n"
        << "Usage:\n"
        << "  " << progName << " genome.fasta reads.fasta [OPTIONS]\n"
        << "\n"
        << "Required arguments:\n"
        << "  genome.fasta            FASTA file containing the reference genome\n"
        << "  reads.fasta             FASTA file containing the reads to align\n"
        << "\n"
        << "Options:\n"
        << "  --kmer_size  INT        Length of k-mers used to index the genome\n"
        << "                            default : 21\n"
        << "                            range   : 1 to 31 (must fit in uint64_t)\n"
        << "  --table_size INT        Number of buckets in the hash table\n"
        << "                            default : 4000000\n"
        << "                            tip     : larger = fewer collisions, more memory\n"
        << "  --heap_size  INT        Maximum number of candidate alignments kept\n"
        << "                          per read in the BestResultsHeap\n"
        << "                            default : 15\n"
        << "  --help                  Print this help message and exit\n"
        << "\n"
        << "Output (TSV on stdout):\n"
        << "  #read  chrom  start  end  mismatches  matches  candidates\n"
        << "  Unmapped reads are reported with chrom=* and numeric fields set to -1.\n"
        << "\n"
        << "Examples:\n"
        << "  " << progName << " genome.fasta reads.fasta\n"
        << "  " << progName << " genome.fasta reads.fasta --kmer_size 15\n"
        << "  " << progName << " genome.fasta reads.fasta --kmer_size 15 --table_size 8000000\n"
        << "  " << progName << " genome.fasta reads.fasta --heap_size 32\n"
        << "\n";
}

int main(int argc, char* argv[]) {
    std::string genomePath;
    std::string readsPath;
    int      k         = 21;
    uint32_t tableSize = 4000000u;
    int      heapSize  = 15;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--help" || arg == "-h") {
            printHelp(argv[0]);
            return 0;
        } else if (arg == "--kmer_size") {
            if (i + 1 >= argc) {
                std::cerr << "Error: --kmer_size requires a value\n";
                return 1;
            }
            k = std::atoi(argv[++i]);
        } else if (arg == "--table_size") {
            if (i + 1 >= argc) {
                std::cerr << "Error: --table_size requires a value\n";
                return 1;
            }
            tableSize = static_cast<uint32_t>(std::atoi(argv[++i]));
        } else if (arg == "--heap_size") {
            if (i + 1 >= argc) {
                std::cerr << "Error: --heap_size requires a value\n";
                return 1;
            }
            heapSize = std::atoi(argv[++i]);
        } else if (arg.size() > 2 && arg[0] == '-' && arg[1] == '-') {
            std::cerr << "Error: unknown option '" << arg << "'\n";
            std::cerr << "Run with --help for usage.\n";
            return 1;
        } else {
            if (genomePath.empty())       genomePath = arg;
            else if (readsPath.empty())   readsPath  = arg;
            else {
                std::cerr << "Error: unexpected argument '" << arg << "'\n";
                std::cerr << "Run with --help for usage.\n";
                return 1;
            }
        }
    }

    if (genomePath.empty() || readsPath.empty()) {
        std::cerr << "Error: genome.fasta and reads.fasta are required.\n";
        std::cerr << "Run with --help for usage.\n";
        return 1;
    }
    if (k < 1 || k > 31) {
        std::cerr << "Error: --kmer_size must be between 1 and 31\n";
        return 1;
    }
    if (tableSize < 1u) {
        std::cerr << "Error: --table_size must be at least 1\n";
        return 1;
    }
    if (heapSize < 1) {
        std::cerr << "Error: --heap_size must be at least 1\n";
        return 1;
    }

    std::cerr << "[load] genome: " << genomePath << "\n";
    GenomeMap genome = readFasta(genomePath);
    std::cerr << "[load] " << genome.size() << " sequence(s)\n";

    HashTable ht(tableSize);
    buildIndex(ht, genome, k);

    std::cerr << "[load] reads: " << readsPath << "\n";
    GenomeMap reads = readFasta(readsPath);
    std::cerr << "[load] " << reads.size() << " read(s)\n";

    mapReads(ht, genome, reads, k, heapSize);
    return 0;
}
