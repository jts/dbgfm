//-----------------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@gmail.com)
// Released under the GPL 
//-----------------------------------------------------
//
// utility - small utility functions that don't fit
// cleanly elsewhere
//
#ifndef UTILITY_H
#define UTILITY_H

#include <string>

// return true if x is a power of 2
#define IS_POWER_OF_2(x) ((x) & ((x) - 1)) == 0

// return the x % y given that y is a power of 2
#define MOD_POWER_2(x, y) (x) & ((y) - 1)

// Calculates the number of bit shifts required to divide a number by divisor
int calculateShiftValue(int divisor);

// convert an integer into its binary encoding 
std::string int2Binary(size_t v, int numBits);

// Complement a base
inline char complement(char base)
{
    switch(base)
    {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: assert(false && "Unknown base!"); return 'N';
    }
}

// reverse-complement a string
inline std::string reverseComplement(const std::string& s)
{
    std::string out(s.length(), 'A');
    size_t last_pos = s.length() - 1;
    for(int i = last_pos; i >= 0; --i)
        out[last_pos - i] = complement(s[i]);
    return out;
}

#endif
