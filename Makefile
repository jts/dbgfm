BAMTOOLS_PATH=/u/jsimpson/simpsonlab/users/jsimpson/software/bamtools
CPP_SRC=fm_index.cpp fm_index_builder.cpp utility.cpp sga_bwt_reader.cpp alphabet.cpp main.cpp
CPP_HEADERS=alphabet.h fm_index_builder.h fm_index.h fm_markers.h \
            huffman_tree_codec.h packed_table_decoder.h sga_bwt_reader.h \
            sga_rlunit.h stream_encoding.h utility.h

# Build the C object files with gcc
%.o: %.c
	gcc -c -o $@ $<

# Build and link the main program
dbgfm: $(CPP_SRC) $(CPP_HEADERS)
		g++ -g -o $@ $(CPP_SRC)

clean:
		rm *.o dbgfm
