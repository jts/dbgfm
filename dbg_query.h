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

namespace DBGQuery
{
    // Returns true if the k-mer string s is a vertex in the
    // de Bruijn graph represented by the provided FM-index
    bool isVertex(const FMIndex* index, const std::string& s);

    // Returns true if the k-mer s has a suffix/prefix extension with the symbol b
    // These functions are reverse-complement aware and will check both strands.
    bool isSuffixNeighbor(const FMIndex* index, const std::string& s, char b);
    bool isPrefixNeighbor(const FMIndex* index, const std::string& s, char b);

    // Returns a string representing the neighboring
    // vertices of the given k-mer s. The neighbors
    // are encoded by the neighbor's extending
    // base in the returned string. For example if the 
    // graph contains these edges:
    //  "ACGT" <-> "CGTA", "ACGT" <-> "CGTG"
    // then getSuffixNeighbors("ACGT") would return "AG".
    // These functions are reverse-complement aware and will check both strands.
    std::string getSuffixNeighbors(const FMIndex* index, const std::string& s);
    std::string getPrefixNeighbors(const FMIndex* index, const std::string& s);
}

#endif
