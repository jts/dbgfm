//-----------------------------------------------
// Copyright 2013 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// DBGQuery - API for querying properties of a
// de Bruijn graph encoded as an FM-index
//
#ifndef DBG_QUERY_H
#define DBG_QUERY_H

#include "fm_index.h"
#include <string>
#include <utility>

namespace DBGQuery
{
    // Returns true if the k-mer string s is a vertex in the
    // de Bruijn graph represented by the provided FM-index
    bool isVertex(const FMIndex* index, const std::string& s);

    // Check for a particular neighbor of k-mer s in the de Bruijn graph.
    // This uses the (k-1) overlap definition of a de Bruijn graph.
    //
    // Let s = aX where X is the (k-1) suffix of s.
    // If the k-mer Xb exists in the FM-index then
    // isSuffixNeighbor(index, s, b) will return true.
    //
    // Similarily let s = Ya where Y is the (k-1) prefix of s.
    // If bY exists in index, isPrefixNeighbor(index, s, b) will return true.
    // These functions are reverse-complement aware so it is sufficient
    // for either Xb or reverseComplement(Xb) to exist to return true.
    bool isSuffixNeighbor(const FMIndex* index, const std::string& s, char b);
    bool isPrefixNeighbor(const FMIndex* index, const std::string& s, char b);

    // Returns a string representing the neighbors of k-mer s. The neighbors
    // are encoded by the single-base extension of the neighboring vertex.
    // For example if the graph contains these edges:
    //  "ACGT" <-> "CGTA", "ACGT" <-> "CGTG"
    // then getSuffixNeighbors(index, "ACGT") would return "AG".
    // Like isSuffixNeighbor/isPrefixNeighbor these functions will check both strands.
    std::string getSuffixNeighbors(const FMIndex* index, const std::string& s);
    std::string getPrefixNeighbors(const FMIndex* index, const std::string& s);

    // Extract a substring of the original text by decompressing a portion
    // of the FM-index. Also return the suffix array index of the
    // substring.
    // idx is the position in the BWT of the last symbol of the substring
    // len is the length of the substring to extract
    std::pair<std::string, size_t>
    extractSubstringAndIndex(const FMIndex* index, size_t idx, size_t len);

    // Extract a substring of the original text.
    // See extractSubstringAndIndex.
    static inline
    std::string
    extractSubstring(const FMIndex* index, size_t idx, size_t len)
    {
        return extractSubstringAndIndex(index, idx, len).first;
    }
}

#endif
