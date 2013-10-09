//-----------------------------------------------
// Copyright 2013 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// DBGQuery - API for querying properties of a
// de Bruijn graph encoded as an FM-index
//
#include <stdio.h>
#include "dbg_query.h"

//
bool DBGQuery::isVertex(const FMIndex* index, const std::string& s)
{
    return index->count(s) > 0 || index->count(reverseComplement(s)) > 0;
}

//
bool DBGQuery::isSuffixNeighbor(const FMIndex* index, const std::string& s, char b)
{
    // Make the neighbor string
    std::string t = s.substr(1) + b;
    return isVertex(index, t);
}

//
bool DBGQuery::isPrefixNeighbor(const FMIndex* index, const std::string& s, char b)
{
    // Make the neighbor string
    std::string t = b + s.substr(0, s.size() - 1);
    return isVertex(index, t);
}

//
std::string DBGQuery::getSuffixNeighbors(const FMIndex* index, const std::string& s)
{
    std::string out;
    for(size_t i = 0; i < 4; ++i)
    {
        char b = "ACGT"[i];
        if(isSuffixNeighbor(index, s, b))
            out.append(1, b);
    }
    return out;
}

//
std::string DBGQuery::getPrefixNeighbors(const FMIndex* index, const std::string& s)
{
    std::string out;
    for(size_t i = 0; i < 4; ++i)
    {
        char b = "ACGT"[i];
        if(isPrefixNeighbor(index, s, b))
            out.append(1, b);
    }
    return out;
}

//
std::pair<std::string, size_t>
DBGQuery::extractSubstringAndIndex(
        const FMIndex* index, size_t idx, size_t len)
{
    std::string out;
    out.reserve(len);

    while(out.length() < len)
    {
        char b = index->getChar(idx);
        if(b == EOF)
            break;

        out.push_back(b);
        idx = index->LF(idx);
    }

    std::reverse(out.begin(), out.end());
    return std::make_pair(out, idx);
}
