//-----------------------------------------------------
// Copyright 2013 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@gmail.com)
// Released under the GPL 
//-----------------------------------------------------
//
// utility - small utility functions that don't fit
// cleanly elsewhere
//
#include <assert.h>
#include "utility.h"

//
int calculateShiftValue(int divisor)
{
    assert(divisor > 0);
    assert(IS_POWER_OF_2(divisor));

    // m_sampleRate is a power of 2, count what bit is set
    unsigned int v = divisor;
    unsigned int c = 0; // c accumulates the total bits set in v

    while(v != 1)
    {
        v >>= 1;
        ++c;
    }
    assert(1 << c == divisor);
    return c;
}

//
std::string int2Binary(size_t v, int numBits)
{
    std::string tmp;
    int bits = sizeof(v) * 8;
    for(int i = bits - 1; i >= 0; --i)
    {
        // test if the i-th bit is set
        size_t mask = 1;
        mask <<= i;
        char b = (v & mask) ? '1' : '0';
        tmp.append(1, b);
    }

    size_t pos = 0;
    if(numBits == 0)
    {
        // Truncate leading zeros
        size_t pos = tmp.find_first_of('1');
        if(pos == std::string::npos)
            pos = tmp.size() - 1;
    }
    else
    {
        pos = tmp.size() - numBits;
    }
    return tmp.substr(pos);
}
