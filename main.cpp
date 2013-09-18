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
    
    // Read the input sequence in its joined form 
    std::string test_file = prefix + ".fa.joined";
    std::ifstream in_file(test_file.c_str());

    std::string sequence;
    getline(in_file, sequence);

    // Count the number of times the reference sequence appears in the
    // FM-index. This must be 1 if the index is correctly loaded.
    printf("//\n// Verifying input sequence is represented by the FM-index\n//\n");
    size_t ref_count = index.count(sequence);
    printf("reference count: %zu\n", ref_count);
    assert(ref_count == 1);
    printf("\n");

    // Test that extract string is working
    printf("//\n// Testing extractSubstring()\n//\n");
    size_t extract_len = sequence.size() > 1000 ? 1000 : sequence.size();
    printf("Extracting last %zu symbols of the test sequence\n", extract_len);
    size_t start_idx = 0;
    std::string extracted = DBGQuery::extractSubstring(&index, start_idx, 1000);
    assert(sequence.find(extracted) != std::string::npos);
    printf("Extracted string matches input sequence\n\n");

    // Test the de bruijn query functions using the graph implied by the reference
    printf("//\n// Testing deBruijn queries for known sequences\n//\n");
    size_t stride = 1000;
    size_t k = 31;
    size_t n_checked = 0;
    size_t n_suffix_branch = 0;
    size_t n_prefix_branch = 0;
    for(size_t idx = 0; idx + k < sequence.size(); idx += stride)
    {
        std::string curr = sequence.substr(idx, k);
        std::string next = sequence.substr(idx + 1, k);

        if(curr.find('$') != std::string::npos || 
           next.find('$') != std::string::npos)
            continue;

        if(curr.size() < k || next.size() < k)
            continue;

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

    printf("num vertices checked: %zu\n", n_checked);
    printf("num suffix branches: %zu\n", n_suffix_branch);
    printf("num prefix branches: %zu\n", n_prefix_branch);
    printf("\n");

    printf("//\n// Testing deBruijn queries for random sequences\n//\n");
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
