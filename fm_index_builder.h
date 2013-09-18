//-----------------------------------------------
// Copyright 2013 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// FMIndexBuilder - Construct an FM-Index from
// an SGA BWT file
//
#ifndef FM_INDEX_BUILDER_H
#define FM_INDEX_BUILDER_H

#include <deque>
#include <fstream>
#include "alphabet.h"
#include "fm_markers.h"
#include "huffman_tree_codec.h"
#include "packed_table_decoder.h"

class FMIndexBuilder
{
    public:
        FMIndexBuilder(const std::string& bwt_filename,
                       size_t small_sample_rate,
                       size_t large_sample_rate);
        ~FMIndexBuilder();
        
        // Get the number of bytes in the compressed string
        size_t getNumStringBytes() const { return m_str_bytes; }

        // Get the number of markers
        size_t getNumSmallMarkers() const { return m_num_small_markers_wrote; }
        size_t getNumLargeMarkers() const { return m_num_large_markers_wrote; }

        size_t getNumStrings() const { return m_strings; }
        size_t getNumSymbols() const { return m_str_symbols; }
        AlphaCount64 getSymbolCounts() const { return m_runningAC; }
        
        // the position in the BWT that represents the full-length string
        size_t getEOFPos() const { return m_eof_pos; }

        // a table to map from huffman symbols to bwt symbols
        PackedTableDecoder getDecoder() const { return m_decoder; }

        // Get the filenames that the result is stored in
        std::string getStringFilename() const;
        std::string getSmallMarkerFilename() const;
        std::string getLargeMarkerFilename() const;
 
    private:
        void build(const std::string& filename);

        void buildSegment(HuffmanTreeCodec<char>& encoder, const std::deque<char>& buffer);
        void buildMarkers();
    
        // the decoding table for the huffman tree we constructed
        PackedTableDecoder m_decoder;

        // the frequency of markers
        size_t m_small_sample_rate;
        size_t m_large_sample_rate;

        // the number of bytes written for the encoded string
        size_t m_str_bytes;

        // the number of symbols written out so far
        size_t m_str_symbols;
        
        // the number of strings in the bwt
        size_t m_strings;
        
        // the FM-index needs to know which BWT symbol corresponds
        // to the full-length string. 
        size_t m_eof_pos;

        // A running count of the number of ACGT$ written
        AlphaCount64 m_runningAC;

        // The last large/small marker that were written out
        SmallMarker m_prevSmallMarker;
        LargeMarker m_prevLargeMarker;

        size_t m_num_small_markers_wrote;
        size_t m_num_large_markers_wrote;

        // temporary output files
        std::ofstream* mp_str_tmp;
        std::ofstream* mp_sm_tmp;
        std::ofstream* mp_lm_tmp;
};

#endif
