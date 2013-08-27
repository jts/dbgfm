dbgfm
=====

An FM-index representation of a de Bruijn graph.

The code in this repository is a stand-alone version of the FM-index from SGA (github/jts/sga).
It is licensed under GPLv3.

## Compiling

The code has no dependencies and should build by just running:
	make

## Testing

To test the code is functioning correctly, you can run:
	make test

This will download human chromosome 20, index it with SGA then perform test queries using dbgfm.

## API

A simple API for querying the structure of the de Bruijn graph is provided. See dbg_query.h.
