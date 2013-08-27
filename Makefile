CPP_SRC := $(wildcard *.cpp)
CPP_HEADERS := $(wildcard *.h)

# SGA is used to prepare the test data
# Point this path to your install
SGA=~/work/code/sga/src/build/SGA/sga

# Build the C object files with gcc
%.o: %.c
	gcc -c -o $@ $<

# Build and link the main program
dbgfm: $(CPP_SRC) $(CPP_HEADERS)
		g++ -O3 -o $@ $(CPP_SRC)

#
# Tests
#

# Download test data
chr20.pp.bwt:
		wget ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37/Primary_Assembly/assembled_chromosomes/FASTA/chr20.fa.gz
		$(SGA) preprocess --permute chr20.fa.gz > chr20.pp.fa
		$(SGA) index --no-reverse chr20.pp.fa

# Run test
test: chr20.pp.bwt dbgfm
		./dbgfm chr20.pp

clean:
		rm *.o dbgfm
