//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// StreamEncoding -- Encode a stream of symbols
// using a huffman encoder. 
//
#ifndef STREAMENCODING_H
#define STREAMENCODING_H

#include "packed_table_decoder.h"
#include "utility.h"

//#define DEBUG_ENCODING 1
#define BITS_PER_BYTE 8

namespace StreamEncode
{

    // Decode functors for the generic decoding function
    struct AlphaCountDecode
    {
        AlphaCountDecode(AlphaCount64& target) : m_target(target) {}
        inline void operator()(int rank)
        {
            m_target.addByIdx(rank, 1);
        }
        AlphaCount64& m_target;
    };

    struct StringDecode
    {
        StringDecode(std::string& target) : m_target(target) {}
        inline void operator()(int rank)
        {
            m_target.append(1, BWT_ALPHABET::getChar(rank));
        }
        std::string& m_target;
    };

    struct BaseCountDecode
    {
        BaseCountDecode(char targetBase, size_t& targetCount) : m_targetBase(targetBase), 
                                                                    m_targetCount(targetCount) {}
        inline void operator()(int rank)
        {
            m_targetCount += (BWT_ALPHABET::getChar(rank) == m_targetBase);
        }
        char m_targetBase;
        size_t& m_targetCount;
    };   

    // Decoder which returns the last base added. This is used to extract a particular character from the stream
    struct SingleBaseDecode
    {
        SingleBaseDecode(char& base) : m_base(base) {}
        inline void operator()(int rank, int /*rl*/)
        {
            m_base = BWT_ALPHABET::getChar(rank);
        }
        char& m_base;
    };

    //
    inline void printEncoding(const std::vector<uint8_t>& output)
    {
        std::cout << "Encoding: ";
        for(size_t i = 0; i < output.size(); ++i)
        {
            std::cout << int2Binary(output[i], 8) << " ";
        }
        std::cout << "\n";
    }

    // Write the code into the output stream starting at currBit
    inline size_t _writeCode(EncodePair& ep, size_t currBit, std::vector<uint8_t>& output)
    {
#ifdef DEBUG_ENCODING
        printEncoding(output);
        std::cout << "Writing the code " << int2Binary(ep.code, ep.bits) << " at bit " << currBit << "\n";
#endif

        size_t code = ep.code;
        int codeBits = ep.bits;
        int bitsRemaining = codeBits;
        int bitOffset = 0;

        while(bitsRemaining > 0)
        {
            // Calculate position to start the write
            int byte = currBit / BITS_PER_BYTE;
            int bitIdx = MOD_POWER_2(currBit, BITS_PER_BYTE);
            int bitsToWrite = std::min((BITS_PER_BYTE - bitIdx), bitsRemaining);

            // Calculate the shift values and masks to apply
            int currPos = (BITS_PER_BYTE - codeBits + bitOffset);
            
            // Mask off the bits we want
            int mask = ((1 << bitsToWrite) - 1) << (codeBits - (bitsToWrite + bitOffset));
            int inCode = code & mask;

#ifdef DEBUG_ENCODING
            std::cout << "Mask: " << int2Binary(mask, codeBits) << "\n";
            std::cout << "Masked: " << int2Binary(inCode,codeBits) << "\n";
#endif

            // Shift the code into position
            if(currPos < bitIdx)
                inCode >>= (bitIdx - currPos);
            else if(currPos > bitIdx)
                inCode <<= (currPos - bitIdx);

#ifdef DEBUG_ENCODING
            std::cout << "Shifted: " << int2Binary(inCode,8) << "\n";
#endif

            // set the value with an OR
            output[byte] |= inCode;

            bitsRemaining -= bitsToWrite;
            bitOffset += bitsToWrite;
            currBit += bitsToWrite;
        }
        return codeBits;
    }

    // Read maxBits from the array starting at currBits and write the value to outCode
    // Returns the number of bits read
    inline void _readCode(const uint16_t currBit, const uint16_t baseShift, const uint16_t mask, const uint16_t input, uint16_t& outCode)
    {
#ifdef DEBUG_ENCODING
        printEncoding(input);
        std::cout << "Reading " << maxBits << " from array starting at " << currBit << "\n";
        std::cout << "Mask " << int2Binary(mask, maxBits) << "\n";
#endif      
        outCode = (input >> (baseShift - currBit)) & mask;
    }
    
    // Encode a stream of characters
    // Returns the number of bytes written
    inline size_t encode(const std::deque<char>& input, const HuffmanTreeCodec<char>& encoder, std::vector<uint8_t>& output)
    {
        // Require the encoder to emit at most 8-bit codes
        assert(encoder.getMaxBits() <= BITS_PER_BYTE);

        // Perform the encoding
        size_t currBit = 0;
        for(size_t i = 0; i < input.size(); ++i)
        {
            assert(currBit / 8 < output.size());

#ifdef DEBUG_ENCODING
            std::cout << "Encoding: " << input[i] << "\n";
#endif
            EncodePair symEP = encoder.encode(input[i]);
            //printf("wcode: %s\n", int2Binary(symEP.code, symEP.bits).c_str());
            currBit += _writeCode(symEP, currBit, output);
        }

        return (currBit + 1) / 8;
    }

    // Decode a stream into the provided functor

#define DECODE_UNIT uint64_t
#define DECODE_UNIT_BYTES sizeof(DECODE_UNIT)
#define DECODE_UNIT_BITS DECODE_UNIT_BYTES * 8
    // Decompress the data starting at pInput. The read cannot exceed the endpoint given by pEnd. Returns
    // the total number of symbols decoded. The out parameters numBitsDecoded is also set.
    template<typename Functor>
    inline size_t decode(const CharPackedTableDecoder& decoder, 
                         const unsigned char* pInput, 
                         const unsigned char* pEnd, 
                         size_t targetSymbols, 
                         DECODE_UNIT& numBitsDecoded, 
                         Functor& functor)
    {
        const std::vector<PACKED_DECODE_TYPE>* p_decode_table = decoder.getTable();
        DECODE_UNIT read_length = decoder.getCodeReadLength();

        // Prime the decode unit by reading bits from the stream
        DECODE_UNIT numBitsBuffered = 0;
        
        DECODE_UNIT decodeUnit = *pInput++;
        for(size_t i = 0; i < DECODE_UNIT_BYTES - 1; ++i)
        {
            decodeUnit <<= BITS_PER_BYTE;
            decodeUnit |= *pInput++;
        }
        numBitsBuffered = BITS_PER_BYTE * DECODE_UNIT_BYTES;
        
        DECODE_UNIT mask = (1 << read_length) - 1;

        // Read data
        numBitsDecoded = 0;
        size_t numSymbolsDecoded = 0;
        while(1)
        {
            // Read a code from the buffered data
            DECODE_UNIT code = 0;
            code = decodeUnit >> (numBitsBuffered - numBitsDecoded - read_length) & mask;
            // Parse the code
            PACKED_DECODE_TYPE packed_code = (*p_decode_table)[code];
            int symbol_rank = UNPACK_SYMBOL(packed_code);
            numBitsDecoded += UNPACK_BITS(packed_code);
            
            //printf("rcode: %s %d\n", int2Binary(code, read_length).c_str(), UNPACK_BITS(packed_code));
            
            // pass the decoded data to the consumer
            functor(symbol_rank);
            numSymbolsDecoded += 1;
            if(numSymbolsDecoded == targetSymbols)
                return numSymbolsDecoded;
            
            // Update the decode unit
            if(numBitsBuffered - numBitsDecoded < 2 * BITS_PER_BYTE)
            {
                for(size_t i = 2; i < DECODE_UNIT_BYTES; ++i)
                    decodeUnit = (decodeUnit << BITS_PER_BYTE) | (pInput <= pEnd ? *pInput++ : 0);
                numBitsBuffered += (BITS_PER_BYTE * (DECODE_UNIT_BYTES - 2));
            }
        }
        return targetSymbols;
    }    
};

#endif
