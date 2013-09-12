//-----------------------------------------------
// Copyright 2013 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// BWTDiskReader - read a bwtdisk file
//
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "bwtdisk_reader.h"

//
BWTDiskReader::BWTDiskReader(const std::string& filename) : m_stage(IOS_HEADER)
{
    m_pReader = new std::ifstream(filename.c_str());
    if(!m_pReader->is_open())
    {
        std::cerr << "Error: could not open " << filename << " for read\n";
        exit(EXIT_FAILURE);
    }
}

//
BWTDiskReader::~BWTDiskReader()
{
    delete m_pReader;
}

//
void BWTDiskReader::discardHeader()
{
    size_t size = 0;
    
    m_pReader->read((char*)&size, sizeof(size));
    m_pReader->read((char*)&m_eof_pos, sizeof(m_eof_pos));

    m_stage = IOS_BWSTR;
    m_num_read = 0;
}

// Read a single base from the BWStr
// The BWT is stored as runs on disk, so this class keeps
// an internal buffer of a single run and emits characters from this buffer
// and performs reads as necessary. If all the runs have been read, emit
// a newline character to signal the end of the BWT
char BWTDiskReader::readChar()
{
    assert(m_stage == IOS_BWSTR);
    
    // Extract a single character from the stream
    char c = m_pReader->get();

    // bwtdisk writes the sentinel position in the header
    // of the file and emits an arbitrary symbol at this location
    // We catch this case here and emit a $ which is what SGA expects.
    if(m_num_read == m_eof_pos)
        c = '$';
    else if(c == EOF) // SGA expects \n at the file's end
        c = '\n';

    m_num_read++;
    return c;
}
