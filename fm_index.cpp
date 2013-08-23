//-----------------------------------------------
// Copyright 2009 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL
//-----------------------------------------------
//
// FMIndex - Run-length encoded Burrows Wheeler transform
//
#include <istream>
#include <queue>
#include <inttypes.h>
#include <stdio.h>
#include "fm_index.h"
#include "utility.h"
#include "sga_bwt_reader.h"
#include "huffman_tree_codec.h"
#include "fm_index_builder.h"

// Parse a BWT from a file
FMIndex::FMIndex(const std::string& filename, int sampleRate) : m_numStrings(0), 
                                                                m_numSymbols(0), 
                                                                m_largeSampleRate(DEFAULT_SAMPLE_RATE_LARGE),
                                                                m_smallSampleRate(sampleRate)
{
    std::cout << "Loading " << filename << "\n";
    loadSGABWT(filename);
    
    /*
    assert(false && "initialize huffman");
    m_rlHuffman = HuffmanTreeCodec<int>();
    m_rlDecodeTable.initialize(m_rlHuffman);

    assert(false && "reader removed");
    IBWTReader* pReader = BWTReader::createReader(filename);
    pReader->read(this);

    //initializeFMIndex();
    delete pReader;
    */
}

void FMIndex::loadSGABWT(const std::string& filename)
{
    FMIndexBuilder builder(filename, m_smallSampleRate, m_largeSampleRate);

}

void FMIndex::setSampleRates(size_t largeSampleRate, size_t smallSampleRate)
{
    m_smallSampleRate = smallSampleRate;
    m_largeSampleRate = largeSampleRate;

    m_smallShiftValue = calculateShiftValue(m_smallSampleRate);
    m_largeShiftValue = calculateShiftValue(m_largeSampleRate);
}

// get the number of markers required to cover the n symbols at sample rate of d
size_t FMIndex::getNumRequiredMarkers(size_t n, size_t d) const
{
    // we place a marker at the beginning (with no accumulated counts), every m_sampleRate
    // bases and one at the very end (with the total counts)
    size_t num_markers = (n % d == 0) ? (n / d) + 1 : (n / d) + 2;
    return num_markers;
}

// Print the BWT
void FMIndex::print() const
{
    assert(false);
}

// Print information about the BWT
void FMIndex::printInfo() const
{
    size_t small_m_size = m_smallMarkers.size() * sizeof(SmallMarker);
    size_t large_m_size = m_largeMarkers.size() * sizeof(LargeMarker);
    size_t total_marker_size = small_m_size + large_m_size;

    size_t bwStr_size = m_string.size();
    size_t other_size = sizeof(*this);
    size_t total_size = total_marker_size + bwStr_size + other_size;

    double mb = (double)(1024 * 1024);
    double total_mb = total_size / mb;
    
    printf("\nFMIndex info:\n");
    printf("Large Sample rate: %zu\n", m_largeSampleRate);
    printf("Small Sample rate: %zu\n", m_smallSampleRate);
    printf("Contains %zu symbols in %zu bytes (%1.4lf symbols per bytes)\n", m_numSymbols, m_string.size(), (double)m_numSymbols / m_string.size());
    printf("Marker Memory -- Small Markers: %zu (%.1lf MB) Large Markers: %zu (%.1lf MB)\n", small_m_size, small_m_size / mb, large_m_size, large_m_size / mb);
    printf("Total Memory -- Markers: %zu (%.1lf MB) Str: %zu (%.1lf MB) Misc: %zu Total: %zu (%lf MB)\n", total_marker_size, total_marker_size / mb, bwStr_size, bwStr_size / mb, other_size, total_size, total_mb);
    printf("N: %zu Bytes per symbol: %lf\n\n", m_numSymbols, (double)total_size / m_numSymbols);
}
