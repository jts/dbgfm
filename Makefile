# Programs
SGA=sga

# Options
CXXFLAGS=-g -O3

# Directories
prefix=/usr/local
bindir=$(prefix)/bin

# Programs to build
bin_PROGRAMS=dbgfm bwtdisk-prepare

# Targets

all: $(bin_PROGRAMS)

clean:
	rm -f $(bin_PROGRAMS) *.o

install: $(bin_PROGRAMS)
	install $(bin_PROGRAMS) $(DESTDIR)$(bindir)

test: $(bin_PROGRAMS) chr20.pp.dbgfm

uninstall:
	cd $(DESTDIR)$(bindir) && rm -f $(bin_PROGRAMS)

.PHONY: all clean install test uninstall
.DELETE_ON_ERROR:
.SECONDARY:

# Build dbgfm

dbgfm_SOURCES = alphabet.cpp bwtdisk_reader.cpp dbg_query.cpp \
	fm_index.cpp fm_index_builder.cpp main.cpp sga_bwt_reader.cpp \
	utility.cpp

dbgfm_HEADERS = alphabet.h bwtdisk_reader.h dbg_query.h fm_index.h \
	fm_index_builder.h fm_markers.h huffman_tree_codec.h \
	packed_table_decoder.h sga_bwt_reader.h sga_rlunit.h \
	stream_encoding.h utility.h

dbgfm: $(dbgfm_SOURCES) $(dbgfm_HEADERS)
	$(CXX) $(INCLUDES) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $(dbgfm_SOURCES) $(LIBS)

# Build bwtdisk-prepare

bwtdisk-prepare: bwtdisk_prepare.cpp
	$(CXX) $(INCLUDES) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

# Tests

chr20.fa.gz:
	wget ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37/Primary_Assembly/assembled_chromosomes/FASTA/chr20.fa.gz

%.pp.fa: %.fa.gz
	$(SGA) preprocess --permute $< >$@

%.bwtdisk: %.fa bwtdisk-prepare
	./run_bwtdisk.sh $<

%.dbgfm: %.bwtdisk dbgfm
	./dbgfm $*
