#!/bin/bash
set -eu

# Prepare the file by concatenating the contigs/sequences into one long string
# separated by $
./bwtdisk-prepare $1 > $1.joined

# Reverse the text, because bwtdisk constructs the BWT of the reverse text by default.
text_rev -vvv $1.joined $1.joined.rev

# Generate the BWT.
bwte -vvv $1.joined.rev

# Rename the bwt
mv $1.joined.rev.bwt `basename $1 .fa`.bwtdisk
