//-----------------------------------------------
// Copyright 2013 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL 
//-----------------------------------------------
//
// bwtdisk_prepare - prepare a FASTA file for indexing
// with bwtdisk
//
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <time.h>

#include "alphabet.h"

int main(int argc, char** argv)
{

    if(argc != 2) {
        fprintf(stderr, "Error: a filename must be provided\n");
        fprintf(stderr, "usage: bwtdisk-prepare filename > output\n");
        exit(EXIT_FAILURE);
    }
    
    std::string filename = argv[1];

    std::string line;
    std::ifstream reader(filename.c_str());
    if(!reader.good()) {
        fprintf(stderr, "Error: could not read %s\n", filename.c_str());
        exit(EXIT_FAILURE);
    }

    srand(time(NULL));

    // Read the fasta file line by line.
    // When we hit a header we output a symbol separating the current record
    // from the last. Non-ACGT symbols in the records cause an error.
    size_t n_records = 0;
    while(getline(reader, line)) {
        if(line.empty())
            continue;
        
        if(line[0] == '>') {

            if(n_records++ > 0)
                std::cout << '$';
        } else {
            std::size_t pos = line.find_first_not_of("ACGT");
            // loop while there are no more ambiguous bases found in line
            while(pos != std::string::npos) {
                if(!IUPAC::isValid(line[pos])) {
                    fprintf(stderr, "Error: invalid IUPAC base found.\n");
                    exit(EXIT_FAILURE);
                }
                std::string possibleSymbols = IUPAC::getPossibleSymbols(line[pos]);
                // choose a random base from the possible symbols
                char base = possibleSymbols[rand() % possibleSymbols.length()];
                line[pos] = base;
                // get the next ambiguous base
                pos = line.find_first_not_of("ACGT");
            }
               
            std::cout << line;
        }
    }

    // Print a final sentinel for the last string
    std::cout << '$';
}
