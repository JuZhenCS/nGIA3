# nGIA3

nGIA3 is a biological-sequence clustering tool. The current implementation
provides the makedb database-packing stage.

## Build and run the C++ packer
The packer requires a C++17 compiler with OpenMP; CUDA is not required.


    cd src
    make
    ./ngia3 makedb -f ../data/go.fas -p go.packed

Run the Python unit tests and compiled C++ parity tests together with:

    make test


The input must end in a newline. Both implementations accept LF or CRLF and
normalize stored records to LF. Records with sequence length greater than or
equal to 65535 are skipped. Retained records are stably sorted by decreasing
sequence length before they are stored.

## Independent Python reference packer

tools/makedb.py independently implements the same 64-bit little-endian packed
format. It uses NumPy automatically when installed and otherwise falls back to
a slower, standard-library-only reference implementation.

    python tools/makedb.py -f data/go.fas -p go.python.packed
    python -m unittest discover -s tests -v

For a supported FASTA input, the C++ and Python output files should have the
same MD5 digest. The Python writer first creates and synchronizes a temporary
file in the destination directory, then atomically replaces the destination.

## Packed file layout

All integer fields are little-endian. Offsets are 64-bit values.

1. One 32-bit retained sequence count.
2. One 32-bit FASTA-name length per retained sequence.
3. One 32-bit sequence length per retained sequence.
4. One packed-data offset per retained sequence.
5. One normalized-FASTA offset per retained sequence.
6. 128 32-bit MinHash signatures per retained sequence.
7. Four bit-plane words per block of 32 residues, prefixed by sequence length.
8. Normalized FASTA records, each using one sequence line and LF delimiters.
