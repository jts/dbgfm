CPP_SRC := $(wildcard *.cpp)
CPP_HEADERS := $(wildcard *.h)

# Build the C object files with gcc
%.o: %.c
	gcc -c -o $@ $<

# Build and link the main program
dbgfm: $(CPP_SRC) $(CPP_HEADERS)
		g++ -g -o $@ $(CPP_SRC)

clean:
		rm *.o dbgfm
