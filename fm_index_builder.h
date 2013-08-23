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

class FMIndexBuilder
{
    public:
        FMIndexBuilder(const std::string& bwt_filename,
                       size_t small_sample_rate,
                       size_t large_sample_rate);
        ~FMIndexBuilder();
 
    private:
        void build(const std::string& filename);

        void buildSegment(HuffmanTreeCodec<char>& encoder, const std::deque<char>& buffer);
        void buildMarkers();

        // the frequency of markers
        size_t m_small_sample_rate;
        size_t m_large_sample_rate;

        // the number of bytes written for the encoded string
        size_t m_str_bytes;

        // the number of symbols written out so far
        size_t m_str_symbols;

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
