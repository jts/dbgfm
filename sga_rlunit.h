//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// RLUnit - A run-length encoded unit of the FM-index
//
#ifndef RLUNIT_H
#define RLUNIT_H

#include <string.h>
#include <inttypes.h>
#include <assert.h>
#include "alphabet.h"

//
#define RL_COUNT_MASK 0x1F  //00011111
#define RL_SYMBOL_MASK 0xE0 //11100000
#define RL_FULL_COUNT 31
#define RL_SYMBOL_SHIFT 5
#define RLE_VALIDATE

// A 5+3 encoded symbol, run pair
// The high 3 bits encodes the symbol to store
// The low 5 bits encodes the length of the run
struct RLUnit
{
    RLUnit() : data(0) {}
    RLUnit(char b) : data(1)
    {
        setChar(b);   
    }

    // Returns true if the count cannot be incremented
    inline bool isFull() const
    {
        return (data & RL_COUNT_MASK) == RL_FULL_COUNT;
    }

    inline bool isEmpty() const
    {
        return (data & RL_COUNT_MASK) == 0;
    }

    inline bool isInitialized() const
    {
        return data > 0;
    }

    // 
    inline void incrementCount()
    {
#ifdef RLE_VALIDATE
        assert(!isFull());
#endif
        ++data;
    }

    // 
    inline void decrementCount()
    {
#ifdef RLE_VALIDATE
        assert(!isEmpty());
#endif
        --data;
    }    

    inline uint8_t getCount() const
    {
#ifdef RLE_VALIDATE
        assert((data & RL_COUNT_MASK) != 0);
#endif
        return data & RL_COUNT_MASK;
    }

    // Set the symbol
    inline void setChar(char symbol)
    {
        // Clear the current symbol
        data &= RL_COUNT_MASK;
        
        uint8_t code = BWT_ALPHABET::getRank(symbol);
        code <<= RL_SYMBOL_SHIFT;
        data |= code;
    }

    // Get the symbol
    inline char getChar() const
    {
        uint8_t code = data & RL_SYMBOL_MASK;
        code >>= RL_SYMBOL_SHIFT;
        return BWT_ALPHABET::getChar(code);
    }

    // 
    uint8_t data;

};

#endif
