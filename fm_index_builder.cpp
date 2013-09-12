//-----------------------------------------------
// Copyright 2013 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// FMIndexBuilder - Construct an FM-Index from
// an SGA BWT file
//
#include <map>
#include "fm_index_builder.h"
#include "sga_bwt_reader.h"
#include "bwtdisk_reader.h"
#include "stream_encoding.h"

FMIndexBuilder::FMIndexBuilder(const std::string& filename, 
                               size_t small_sample_rate,
                               size_t large_sample_rate)
{
    // Create temporary files for the 3 components of the index
    mp_str_tmp = new std::ofstream(getStringFilename().c_str(), std::ios::binary);
    mp_sm_tmp = new std::ofstream(getSmallMarkerFilename().c_str(), std::ios::binary);
    mp_lm_tmp = new std::ofstream(getLargeMarkerFilename().c_str(), std::ios::binary);

    m_small_sample_rate = small_sample_rate;
    m_large_sample_rate = large_sample_rate;


    build(filename);
}

FMIndexBuilder::~FMIndexBuilder()
{
    assert(mp_str_tmp == NULL);
    assert(mp_sm_tmp == NULL);
    assert(mp_lm_tmp == NULL);
}

std::string FMIndexBuilder::getStringFilename() const
{
    return "dbgfm.str";
}

std::string FMIndexBuilder::getSmallMarkerFilename() const
{
    return "dbgfm.sm";

}

std::string FMIndexBuilder::getLargeMarkerFilename() const
{
    return "dbgfm.lm";
}

void FMIndexBuilder::build(const std::string& filename)
{
    // Initialization
    m_str_bytes = 0;
    m_str_symbols = 0;
    m_num_large_markers_wrote = 0;
    m_num_small_markers_wrote = 0;

    //
    // Step 1: make a symbol -> count map and use it to build a huffman tree
    //
    std::map<char, size_t> count_map;
    BWTDiskReader* p_reader = new BWTDiskReader(filename);

    // Discard header for now
    p_reader->discardHeader();

    // Read one symbol from the bwt at a time
    char b;

    while((b = p_reader->readChar()) != '\n') {
        count_map[b]++;
    }

    HuffmanTreeCodec<char> encoder(count_map);
    m_decoder.initialize(encoder);
    m_strings = count_map['$'];

    for(std::map<char, size_t>::iterator iter = count_map.begin();
        iter != count_map.end(); ++iter) {
        printf("%c %zu\n", iter->first, iter->second);
    }

    std::cout << "Bits required for string: " << encoder.getRequiredBits(count_map) << "\n";

    //
    // Step 2: use the huffman tree to compress the string
    //

    // re-initialize the reader
    delete p_reader;
    p_reader = new BWTDiskReader(filename);
    p_reader->discardHeader();

    // We buffer 128 or 256 symbols at a time and huffman-encode each segment
    std::deque<char> buffer;

    while((b = p_reader->readChar()) != '\n')
    {
        buffer.push_back(b);
        if(buffer.size() == m_small_sample_rate)
        {
            buildSegment(encoder, buffer);
            buffer.clear();
        }
    }
    
    // Build a segment for the remaining symbols
    if(!buffer.empty())
        buildSegment(encoder, buffer);

    delete p_reader;

    delete mp_str_tmp;
    delete mp_sm_tmp;
    delete mp_lm_tmp;

    mp_str_tmp = NULL;
    mp_sm_tmp = NULL;
    mp_lm_tmp = NULL;
}

void FMIndexBuilder::buildSegment(HuffmanTreeCodec<char>& encoder,
                                  const std::deque<char>& buffer)
{
    // output any markers needed
    buildMarkers();
    
    // Update the occurrence counts for the incoming symbols
    for(size_t i = 0; i < buffer.size(); ++i)
        m_runningAC.increment(buffer[i]);

    // make a buffer that is large enough to store the encoded data in the worst case
    size_t max_bits = encoder.getMaxBits() * buffer.size();
    size_t max_bytes = max_bits / 8;
    std::vector<uint8_t> output(max_bytes, 0);

    size_t bytes = StreamEncode::encode(buffer, encoder, output);
    mp_str_tmp->write(reinterpret_cast<const char*>(&output[0]), bytes);

    size_t bits_read = 0;
    PackedTableDecoder decoder;
    decoder.initialize(encoder);
    std::string str;
    StreamEncode::StringDecode sd(str);
    StreamEncode::decode(decoder, &output[0], &output[0] + output.size() - 1, buffer.size(), bits_read, sd);

    std::string e;
    for(size_t i = 0; i < buffer.size(); ++i)
        e.append(1, buffer[i]);

    assert(e == str);
//    printf("E: %s\n", e.c_str());
//    printf("D: %s\n", str.c_str());

    m_str_symbols += buffer.size();
    m_str_bytes += bytes;
}

void FMIndexBuilder::buildMarkers()
{
    size_t starting_byte = m_str_bytes;

    // Do we need to place new large markers?
    while((m_str_symbols / m_large_sample_rate) + 1 > m_num_large_markers_wrote)
    {
        // Build a new large marker with the accumulated counts up to this point
        LargeMarker marker;
        marker.byteIndex = starting_byte;
        marker.counts = m_runningAC;
        m_prevLargeMarker = marker;

        // Write the marker to the temp file
        mp_lm_tmp->write(reinterpret_cast<const char*>(&m_prevLargeMarker), sizeof(LargeMarker));
        m_num_large_markers_wrote += 1;
    }

    // We place a new SmallMarkers for every segment. 
    AlphaCount16 smallAC;
    for(size_t j = 0; j < BWT_ALPHABET::size; ++j)
    {
        size_t v = m_runningAC.getByIdx(j) - m_prevLargeMarker.counts.getByIdx(j);
        if(v > smallAC.getMaxValue())
        {
            std::cerr << "Error: Number of symbols exceeds the maximum value (" << v << " > " << smallAC.getMaxValue() << ")\n";
            std::cerr << "RunningAC: " << m_runningAC << "\n";
            std::cerr << "PrevAC: " << m_prevLargeMarker.counts << "\n";
            std::cerr << "SmallAC:" << smallAC << "\n";
            exit(EXIT_FAILURE);
        }
        smallAC.setByIdx(j, v);
    }
    
    // Construct the small marker
    SmallMarker smallMarker;
    smallMarker.byteCount = starting_byte - m_prevLargeMarker.byteIndex;
    smallMarker.counts = smallAC;        
    m_prevSmallMarker = smallMarker;

    // write it to disk
    mp_sm_tmp->write(reinterpret_cast<const char*>(&m_prevSmallMarker), sizeof(SmallMarker));
    m_num_small_markers_wrote += 1;
}
