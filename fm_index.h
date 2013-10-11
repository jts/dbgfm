//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// FMIndex - Run-length encoded Burrows Wheeler transform
//
#ifndef FMIndex_H
#define FMIndex_H

#include <deque>
#include <utility>
#include "fm_markers.h"
#include "stream_encoding.h"
#include "packed_table_decoder.h"

// Defines
#define FMINDEX_VALIDATE 1

typedef std::vector<uint8_t> FMBytes;

//
// FMIndex
//
class FMIndex
{
    public:
    
        // Constructors
        FMIndex(const std::string& filename, int sampleRate = DEFAULT_SAMPLE_RATE_SMALL);

        // test that the FM-index is correctly initialized
        // by checking against the on-disk bwt
        void verify(const std::string& bwt_filename);

        //    
        void setSampleRates(size_t largeSampleRate, size_t smallSampleRate);
        void initializeFMIndex(AlphaCount64& running_ac);

        // Update the suffix array interval and return whether the
        // interval is open
        bool updateInterval(size_t& lower, size_t& upper, char c) const
        {
                size_t p = getPC(c);
                lower = p + getOcc(c, lower - 1);
                upper = p + getOcc(c, upper) - 1;
                return lower <= upper;
        }

        // Return the suffix array interval of the string
        std::pair<size_t, size_t> findInterval(const std::string& s) const
        {
            assert(!s.empty());
            int j = s.size() - 1;
            char curr = s[j];
            
            // Initialize interval to the range for the last symbol of s
            size_t lower = getPC(curr);
            size_t upper = lower + getOcc(curr, getBWLen() - 1) - 1;

            --j;
            for(;j >= 0; --j)
            {
                curr = s[j];
                // update interval
                size_t p = getPC(curr);
                lower = p + getOcc(curr, lower - 1);
                upper = p + getOcc(curr, upper) - 1;

                if(lower > upper)
                    return std::make_pair(lower, lower - 1);
            }
            return std::make_pair(lower, upper);
        }

        // Count the number of occurrences of the string s in the original text
        size_t count(const std::string& s) const
        {
            std::pair<size_t, size_t> x = findInterval(s);
            return x.second - x.first + 1;
        }

        // Perform the LF mapping
        // Let SA[idx] = i.
        // This function returns idx'
        // where SA[idx'] = i - 1.
        //
        // This function will assert if idx == m_eof_pos
        inline size_t LF(size_t idx) const
        {
            char b = getChar(idx);
            assert(b != EOF);
            size_t p = getPC(b);
            size_t o = idx > 0 ? getOcc(b, idx - 1) : 0;
            return p + o;
        }
        
        // Returns BWT[idx]. 
        // This function will return EOF when SA[idx] == 0
        inline char getChar(size_t idx) const
        {
            // Decompress stream up to the (idx + 1) character and return the last decompressed symbol
            const LargeMarker marker = getLowerMarker(idx);
            size_t current_position = marker.getActualPosition();
            size_t numToCount = idx - current_position + 1;
            //assert(numToCount < m_smallSampleRate);
            size_t symbol_index = marker.byteIndex;
            DECODE_UNIT numBitsRead = 0;

            char outBase;
            StreamEncode::SingleBaseDecode sbd(outBase);
            StreamEncode::decode(m_decoder, &m_string[symbol_index], &m_string.back(), numToCount, numBitsRead, sbd);
            return idx != m_eof_pos ? outBase : EOF;
        }

        // Get the greatest interpolated marker whose position is less than or equal to position
        inline LargeMarker getLowerMarker(size_t position) const
        {
            size_t target_small_idx = position >> m_smallShiftValue;
            return getInterpolatedMarker(target_small_idx);
        }

        // Return a LargeMarker with values that are interpolated by adding
        // the relative count nearest to the requested position to the last
        // LargeMarker
        inline LargeMarker getInterpolatedMarker(size_t target_small_idx) const
        {
            // Calculate the position of the LargeMarker that the SmallMarker is relative to
            size_t target_position = target_small_idx << m_smallShiftValue;
            size_t curr_large_idx = target_position >> m_largeShiftValue;
            LargeMarker absoluteMarker = m_largeMarkers[curr_large_idx];
            assert(target_small_idx < m_smallMarkers.size());
            const SmallMarker& relative = m_smallMarkers[target_small_idx];
            alphacount_add16(absoluteMarker.counts, relative.counts);
            absoluteMarker.byteIndex += relative.byteCount;
            return absoluteMarker;
        }

        inline size_t getPC(char b) const { return m_predCount.get(b); }

        // Return the number of times char b appears in bwt[0, idx]
        inline size_t getOcc(char b, size_t idx) const
        {
            // The counts in the marker are not inclusive so we increment the index by 1.
            ++idx;

            const LargeMarker marker = getLowerMarker(idx);
            size_t current_position = marker.getActualPosition();
            size_t numToCount = idx - current_position;
            assert(numToCount < m_smallSampleRate);
            size_t running_count = marker.counts.get(b);
            size_t symbol_index = marker.byteIndex;
            StreamEncode::BaseCountDecode bcd(b, running_count);
            DECODE_UNIT numBitsRead = 0;
            StreamEncode::decode(m_decoder, &m_string[symbol_index], &m_string.back(), numToCount, numBitsRead, bcd);
            // The EOF marker symbol is stored in the BWT as a '$'.
            // Subtract one from the count of '$' when the index is
            // larger than the position of the EOF marker.
            if (b == '$' && idx > m_eof_pos)
                --running_count;
            return running_count;
        }

        // Return the number of times each symbol in the alphabet appears in bwt[0, idx]
        inline AlphaCount64 getFullOcc(size_t idx) const 
        { 
            // The counts in the marker are not inclusive so we increment the index by 1.
            ++idx;

            const LargeMarker marker = getLowerMarker(idx);
            size_t current_position = marker.getActualPosition();
            AlphaCount64 running_count = marker.counts;
            size_t numToCount = idx - current_position;

            assert(numToCount < m_smallSampleRate);
            size_t symbol_index = marker.byteIndex;
            StreamEncode::AlphaCountDecode acd(running_count);
            DECODE_UNIT numBitsRead = 0;
            StreamEncode::decode(m_decoder, &m_string[symbol_index], &m_string.back(), numToCount, numBitsRead, acd);
            return running_count;
        }

        // Return the number of times each symbol in the alphabet appears ins bwt[idx0, idx1]
        inline AlphaCount64 getOccDiff(size_t idx0, size_t idx1) const 
        { 
            return getFullOcc(idx1) - getFullOcc(idx0); 
        }

        inline size_t getNumStrings() const { return m_numStrings; } 
        inline size_t getBWLen() const { return m_numSymbols; }
        inline size_t getNumBytes() const { return m_string.size(); }
        inline size_t getSmallSampleRate() const { return m_smallSampleRate; }

        // Return the first letter of the suffix starting at idx
        inline char getF(size_t idx) const
        {
            size_t ci = 0;
            while(ci < BWT_ALPHABET::size && m_predCount.getByIdx(ci) <= idx)
                ci++;
            assert(ci != 0);
            return BWT_ALPHABET::getChar(ci - 1);
        }

        // Print the size of the BWT
        void printInfo() const;
        void print() const;
        void printRunLengths() const;
        
        void decodeToFile(const std::string& file);

        // IO
        friend class BWTReaderBinary;
        friend class BWTWriterBinary;
        friend class BWTReaderAscii;
        friend class BWTWriterAscii;

        // Default sample rates for the large (64-bit) and small (8-bit) occurrence markers
        static const int DEFAULT_SAMPLE_RATE_LARGE = 16384;
        static const int DEFAULT_SAMPLE_RATE_SMALL = 128;

    private:


        // Default constructor is not allowed
        FMIndex() {}
        
        // Load an SGA-encoded bwt
        void loadBWT(const std::string& filename);

        // this class consumes huffman codes and emits the symbols they represent
        PackedTableDecoder m_decoder;

        // The C(a) array
        AlphaCount64 m_predCount;
        
        // The compressed bwt string
        FMBytes m_string;

        // The marker vectors
        LargeMarkerVector m_largeMarkers;
        SmallMarkerVector m_smallMarkers;

        // The number of strings in the collection
        size_t m_numStrings;

        // The total length of the bw string
        size_t m_numSymbols;
        
        // Within our BWT implemention the symbol
        // for the full-length suffix is encoded with a $.
        // When extracting the original text from the BWT, 
        // we store this position so we know when to stop extracting.
        // getChar(idx) will return EOF when idx == m_eof_pos
        size_t m_eof_pos;

        // The sample rate used for the markers
        size_t m_largeSampleRate;
        size_t m_smallSampleRate;

        // The amount to shift values by to divide by m_sampleRate
        int m_smallShiftValue;
        int m_largeShiftValue;

};
#endif
