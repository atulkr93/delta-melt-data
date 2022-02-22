#!/bin/csh

#
# Basic 2D Phase-Sensitive Processing:
#   Cosine-Bells are used in both dimensions.
#   Use of "ZF -auto" doubles size, then rounds to power of 2.
#   Use of "FT -auto" chooses correct Transform mode.
#   Imaginaries are deleted with "-di" in each dimension.
#   Phase corrections should be inserted by hand.

nmrPipe -in test.fid \
#| nmrPipe  -fn SOL \
#| nmrPipe  -fn POLY -time -avg\
#| nmrPipe  -fn GM -g1 5 -g2 15 -g3 0.0 -c 1 \
| nmrPipe  -fn SP -off 0.45 -end 0.98 -pow 1 -c 0.5   \
| nmrPipe  -fn ZF -size 4096                               \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 -12.0 -p1 0.00 -di -verb         \
#| nmrPipe  -fn POLY -auto \
| nmrPipe  -fn EXT -x1 2.9ppm -xn 2.6ppm -sw -verb \
| nmrPipe  -fn TP                                     \
#| nmrPipe  -fn GM -g1 5 -g2 15 -g3 0.0 -c 1 \
| nmrPipe  -fn SP -off 0.45 -end 0.98 -pow 1 -c 1.0    \
| nmrPipe  -fn ZF -size 4096                               \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 -90 -p1 180.00 -di -verb         \
| nmrPipe  -fn EXT -x1 13.2ppm -xn 6.8ppm -sw -verb \
| nmrPipe  -fn TP                                     \
   -ov -out test.ft2
#pipe2ucsf test.ft2 A6DNAm3T9_150ms_NOESY.ucsf
