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
                                                                m_numSymbols(0)
{
    setSampleRates(DEFAULT_SAMPLE_RATE_LARGE, sampleRate);

    std::cout << "Loading " << filename << "\n";
    loadBWT(filename);
}

//
void FMIndex::loadBWT(const std::string& filename)
{
    FMIndexBuilder builder(filename, m_smallSampleRate, m_largeSampleRate);

    size_t n = 0;

    // Load the compressed string from the file
    std::ifstream str_reader(builder.getStringFilename().c_str());
    n = builder.getNumStringBytes();
    m_string.resize(n);
    str_reader.read(reinterpret_cast<char*>(&m_string[0]), n);

    // Load the small markers from the file
    std::ifstream sm_reader(builder.getSmallMarkerFilename().c_str());
    n = builder.getNumSmallMarkers();
    m_smallMarkers.resize(n);
    sm_reader.read(reinterpret_cast<char*>(&m_smallMarkers[0]), sizeof(SmallMarker) * n);
    
    // Load the large markers from the file
    std::ifstream lm_reader(builder.getLargeMarkerFilename().c_str());
    n = builder.getNumLargeMarkers();
    m_largeMarkers.resize(n);
    lm_reader.read(reinterpret_cast<char*>(&m_largeMarkers[0]), sizeof(LargeMarker) * n);

    m_numStrings = builder.getNumStrings();
    m_numSymbols = builder.getNumSymbols();

    AlphaCount64 totals = builder.getSymbolCounts();
    assert(totals.get('$') + 
           totals.get('A') + 
           totals.get('C') + 
           totals.get('G') + 
           totals.get('T') == m_numSymbols);

    m_predCount.set('$', 0);
    m_predCount.set('A', totals.get('$')); 
    m_predCount.set('C', m_predCount.get('A') + totals.get('A'));
    m_predCount.set('G', m_predCount.get('C') + totals.get('C'));
    m_predCount.set('T', m_predCount.get('G') + totals.get('G'));
    assert(m_predCount.get('T') + totals.get('T') == m_numSymbols);

    m_decoder = builder.getDecoder();
    m_eof_pos = builder.getEOFPos();

    printInfo();
}

//
void FMIndex::setSampleRates(size_t largeSampleRate, size_t smallSampleRate)
{
    m_smallSampleRate = smallSampleRate;
    m_largeSampleRate = largeSampleRate;

    m_smallShiftValue = calculateShiftValue(m_smallSampleRate);
    m_largeShiftValue = calculateShiftValue(m_largeSampleRate);
}

// Verify that the index is set up correctly
// by comparing it to the on-disk version.
// This is SLOW
void FMIndex::verify(const std::string& filename)
{
    SGABWTReader* p_reader = new SGABWTReader(filename);

    // Discard header for now
    size_t n1, n2;
    BWFlag flag;
    p_reader->readHeader(n1, n2, flag);


    AlphaCount64 running_count;
    // Read one symbol from the bwt at a time
    size_t i = 0;
    char b;
    while((b = p_reader->readChar()) != '\n')
    {
        // Verify that the symbol at position i matches the symbol
        // read from the disk
        char s = getChar(i);
        assert(s == b);

        // Verifiy that the counts interpolated from the markers
        // are correct
        running_count.increment(b);
        
        // single symbol count
        size_t occ = getOcc(s, i);
        assert(occ == running_count.get(b));
        
        // full count
        AlphaCount64 full_occ = getFullOcc(i);
        assert(full_occ == running_count);

        i++;
    }
    printf("Verified all %zu match expected\n", i);
    delete p_reader;
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
    printf("Contains %zu symbols in %zu bytes (%1.4lf symbols per byte)\n", m_numSymbols, m_string.size(), (double)m_numSymbols / m_string.size());
    printf("Marker Memory -- Small Markers: %zu (%.1lf MB) Large Markers: %zu (%.1lf MB)\n", small_m_size, small_m_size / mb, large_m_size, large_m_size / mb);
    printf("Total Memory -- Markers: %zu (%.1lf MB) Str: %zu (%.1lf MB) Misc: %zu Total: %zu (%lf MB)\n", total_marker_size, total_marker_size / mb, bwStr_size, bwStr_size / mb, other_size, total_size, total_mb);
    printf("N: %zu Bytes per symbol: %lf\n\n", m_numSymbols, (double)total_size / m_numSymbols);
}
