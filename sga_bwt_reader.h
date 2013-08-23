//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// SGABWTReader - read sga's bwt file
//
#ifndef SGABWTREADER_H
#define SGABWTREADER_H

#include <fstream>
#include "sga_rlunit.h"

// State enum to track what segment of the 
// BWT file is being parsed 
enum BWIOStage
{
    IOS_NONE,
    IOS_HEADER,
    IOS_BWSTR,
    IOS_PC,
    IOS_OCC,
    IOS_DONE
};

//
enum BWFlag
{
    BWF_NOFMI = 0,
    BWF_HASFMI
};

// Magic number that the sga file starts with
const uint16_t RLBWT_FILE_MAGIC = 0xCACA;

class SGABWTReader
{
    public:
        SGABWTReader(const std::string& filename);
        ~SGABWTReader();

        void readHeader(size_t& num_strings, size_t& num_symbols, BWFlag& flag);
        char readChar();

    private:
        BWIOStage m_stage;
        std::ifstream* m_pReader;
        RLUnit m_currRun;
        size_t m_numRunsOnDisk;
        size_t m_numRunsRead;
};

#endif
