#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include "fm_index.h"
#include "dbg_query.h"

// Return a random string of length n
std::string getRandomSequence(size_t n)
{
    std::string o;
    for(size_t i = 0; i < n; i++)
        o.append(1, "ACGT"[rand() % 4]);
    return o;
}


int main(int argc, char** argv)
{
    if(argc != 2)
    {
        printf("usage: ./dbgfm <reference_prefix>\n");
        exit(EXIT_FAILURE);
    }

    printf("Loading FM-index\n");
    std::string prefix = argv[1];
    std::string test_bwt = prefix + ".bwtdisk";
    FMIndex index(test_bwt, 256);

    // Verify that the FM-index data structures are set correctly
    //index.verify(test_bwt);

    // Read the fasta sequence as a single string
    std::string test_fa = prefix + ".fa";
    printf("Loading %s\n", test_fa.c_str());
    std::ifstream in_fasta(test_fa.c_str());

    std::string header;
    std::string sequence;
    getline(in_fasta, header);
    getline(in_fasta, sequence);

    // Count the number of times the reference sequence appears in the
    // FM-index. This must be 1 if the index is correctly loaded.
    printf("Verifying reference sequence is represented by FM-index\n");
    size_t ref_count = index.count(sequence);
    printf("\treference count: %zu\n", ref_count);
    assert(ref_count == 1);

    // Test the de bruijn query functions using the graph implied by the reference
    size_t stride = 1000;
    size_t k = 31;
    size_t n_checked = 0;
    size_t n_suffix_branch = 0;
    size_t n_prefix_branch = 0;
    for(size_t idx = 0; idx < sequence.size() - k; idx += stride)
    {
        std::string curr = sequence.substr(idx, k);
        std::string next = sequence.substr(idx + 1, k);

        // Check whether these k-mers are in the vertex set
        assert(DBGQuery::isVertex(&index, curr));
        assert(DBGQuery::isVertex(&index, next));
        
        // Check reverse-complements as well
        assert(DBGQuery::isVertex(&index, reverseComplement(curr)));
        assert(DBGQuery::isVertex(&index, reverseComplement(next)));

        // Check whether there is an edge between the k-mers in the implicit graph
        char c_extend = curr[0];
        char n_extend = next[next.size() - 1];

        assert(DBGQuery::isSuffixNeighbor(&index, curr, n_extend));
        assert(DBGQuery::isPrefixNeighbor(&index, next, c_extend));

        // Check the complete neighbor set for each vertex
        std::string c_neighbors = DBGQuery::getSuffixNeighbors(&index, curr);
        std::string n_neighbors = DBGQuery::getPrefixNeighbors(&index, next);

        assert(c_neighbors.find_first_of(n_extend) != std::string::npos);
        assert(n_neighbors.find_first_of(c_extend) != std::string::npos);

        n_checked += 1;
        n_suffix_branch += c_neighbors.size() > 1;
        n_prefix_branch += n_neighbors.size() > 1;

        if(n_checked % 1000 == 0)
            printf("Checked %zu vertices in the dbg graph [curr idx: %zu]\n", n_checked, idx);
    }

    printf("Done reference graph checks\n");
    printf("\tnum vertices checked: %zu\n", n_checked);
    printf("\tnum suffix branches: %zu\n", n_suffix_branch);
    printf("\tnum prefix branches: %zu\n", n_prefix_branch);

    // Check whether random strings are vertices in the graph
    for(size_t k = 11; k <= 31; k += 5)
    {
        printf("Performing random vertex queries for %zu-mers\n", k);
        size_t n_random_checked = 0;
        size_t n_random_passed = 0;
        for(size_t n = 0; n < 10000; ++n)
        {
            std::string r = getRandomSequence(k);
            n_random_checked += 1;
            n_random_passed += DBGQuery::isVertex(&index, r);
        }

        printf("\tnum checked: %zu\n", n_random_checked);
        printf("\tnum in graph: %zu (%.3lf)\n", n_random_passed, (double)n_random_passed / n_random_checked);
    }
}
