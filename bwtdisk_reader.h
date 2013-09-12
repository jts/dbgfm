//-----------------------------------------------
// Copyright 2013 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// BWTDiskReader - read a bwtdisk file
//
#ifndef BWTDISK_READER_H
#define BWTDISK_READER_H

#include <fstream>
#include "sga_bwt_reader.h"

class BWTDiskReader
{
    public:
        BWTDiskReader(const std::string& filename);
        ~BWTDiskReader();

        void discardHeader();
        char readChar();

    private:
        std::ifstream* m_pReader;
        BWIOStage m_stage;
        size_t m_eof_pos;
        size_t m_num_read;
};

#endif
