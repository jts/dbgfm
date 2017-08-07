PACKAGE=dbgfm

# Programs
SGA=sga

# Options
CXXFLAGS=-g -O3

# Directories
prefix=/usr/local
bindir=$(prefix)/bin
includedir=$(prefix)/include
libdir=$(prefix)/lib
pkgincludedir=$(includedir)/$(PACKAGE)

# Programs and libraries to build
PROGRAMS=dbgfm bwtdisk-prepare
LIBRARIES=libdbgfm.a

# Targets

all: $(PROGRAMS)

clean:
	rm -f $(LIBRARIES) $(PROGRAMS) *.o

install: $(PROGRAMS)
	install $(PROGRAMS) $(DESTDIR)$(bindir)
	install $(LIBRARIES) $(DESTDIR)$(libdir)
	install -d $(DESTDIR)$(pkgincludedir)
	install $(HEADERS) $(DESTDIR)$(pkgincludedir)

test: $(PROGRAMS) chr20.pp.dbgfm

uninstall:
	-cd $(DESTDIR)$(bindir) && rm -f $(PROGRAMS)
	-cd $(DESTDIR)$(libdir) && rm -f $(LIBRARIES)
	-cd $(DESTDIR)$(pkgincludedir) && rm -f $(HEADERS)
	-rmdir $(DESTDIR)$(pkgincludedir)

.PHONY: all clean install test uninstall
.DELETE_ON_ERROR:
.SECONDARY:

# Headers

HEADERS = alphabet.h bwtdisk_reader.h dbg_query.h fm_index.h \
	fm_index_builder.h fm_markers.h huffman_tree_codec.h \
	packed_table_decoder.h sga_bwt_reader.h sga_rlunit.h \
	stream_encoding.h utility.h

# Build libdbgfm.a

libdbgfm_a_OBJECTS = alphabet.o bwtdisk_reader.o dbg_query.o \
	fm_index.o fm_index_builder.o sga_bwt_reader.o utility.o

libdbgfm.a: $(libdbgfm_a_OBJECTS) $(HEADERS)
	$(AR) crs $@ $(libdbgfm_a_OBJECTS)

# Build dbgfm

dbgfm: main.o libdbgfm.a
	$(CXX) $(INCLUDES) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

# Build bwtdisk-prepare

bwtdisk-prepare: bwtdisk_prepare.o
	$(CXX) $(INCLUDES) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS) -L. -ldbgfm 

# Tests

chr20.fa.gz:
	wget ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37/Primary_Assembly/assembled_chromosomes/FASTA/chr20.fa.gz

%.pp.fa: %.fa.gz
	$(SGA) preprocess --permute $< >$@

%.bwtdisk: %.fa bwtdisk-prepare
	./run_bwtdisk.sh $<

%.dbgfm: %.bwtdisk dbgfm
	./dbgfm $*
