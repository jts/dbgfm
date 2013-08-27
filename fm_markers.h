//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// FMMarkers - Marker classes used in the FM-index
// implementation
//
#ifndef FMMARKERS_H
#define FMMARKERS_H

#include <vector>
#include "alphabet.h"

// LargeMarker - To allow random access to the 
// BWT symbols and implement the occurrence array
// we keep a vector of symbol counts every D1 symbols.
// These counts are the absolute number of times each
// symbol has been seen up to that point.
// 
struct LargeMarker
{
    LargeMarker() : byteIndex(0) {}

    // Calculate the actual position in the uncompressed BWT of this marker
    // This is the number of symbols preceding this marker
    inline size_t getActualPosition() const
    {
        return counts.getSum();
    }

    void print() const
    {
        std::cout << "Large marker actual pos: " << getActualPosition() << "\n";
        std::cout << "Marker byte index: " << byteIndex << "\n";
        std::cout << "Marker counts: ";
        for(int i = 0; i < BWT_ALPHABET::size; ++i)
            std::cout << (int)counts.getByIdx(i) << " ";
        std::cout << "\n";
    }    

    // Returns true if the data in the markers is identical
    bool operator==(const LargeMarker& rhs)
    {
        for(size_t i = 0; i < BWT_ALPHABET::size; ++i)
        {
            if(counts.getByIdx(i) != rhs.counts.getByIdx(i))
                return false;
        }
        return byteIndex == rhs.byteIndex;
    }

    // The number of times each symbol has been seen up to this marker
    AlphaCount64 counts; 

    // The index in the compressed string the symbols that starts after
    // this marker. That is, if C = getActualPosition(), then
    // the compressed data containing the B[C] is at byteIndex. This is not necessarily
    // a valid index if there is a marker after the last symbol in the BWT.
    size_t byteIndex;
};
typedef std::vector<LargeMarker> LargeMarkerVector;

// SmallMarker - Small markers contain the counts
// within an individual block of the BWT. In other words
// the small marker contains the count for the last D2 symbols
// 
struct SmallMarker
{
    SmallMarker() : byteCount(0) {}

    // Calculate the actual position in the uncompressed BWT of this marker
    // This is the number of symbols preceding this marker
    inline size_t getCountSum() const
    {
        return counts.getSum();
    }

    void print() const
    {
        for(int i = 0; i < BWT_ALPHABET::size; ++i)
            std::cout << (int)counts.getByIdx(i) << " ";
        std::cout << "\n";
    }

    // The number of times each symbol has been seen up to this marker
    AlphaCount16 counts; 

    // The number of compressed bytes in this block
    uint16_t byteCount;
};
typedef std::vector<SmallMarker> SmallMarkerVector;

#endif
