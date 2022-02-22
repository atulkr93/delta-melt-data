#!/bin/csh

#
# Basic 2D Phase-Sensitive Processing:
#   Cosine-Bells are used in both dimensions.
#   Use of "ZF -auto" doubles size, then rounds to power of 2.
#   Use of "FT -auto" chooses correct Transform mode.
#   Imaginaries are deleted with "-di" in each dimension.
#   Phase corrections should be inserted by hand.

nmrPipe -in test.fid \
| nmrPipe  -fn  POLY -time \
| nmrPipe  -fn  SOL \
| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 2 -c 0.5    \
| nmrPipe  -fn ZF -size 8192                              \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 64.00 -p1 0.00 -di -verb         \
#| nmrPipe  -fn EXT -xn 5.8ppm -x1 6.05ppm -sw -verb    \
| nmrPipe  -fn EXT -xn 5.45ppm -x1 6.15ppm -sw -verb    \
| nmrPipe  -fn TP                                     \
| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 2 -c 0.5    \
| nmrPipe  -fn ZF -size 4096                               \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 120.00 -p1 -87.00 -di -verb         \
| nmrPipe  -fn EXT -xn 83.2ppm -x1 86.2ppm -sw -verb    \
#| nmrPipe  -fn EXT -xn 84ppm -x1 86ppm -sw -verb    \
   -ov -out test.ft2
#pipe2txt.tcl test.ft2 > see.txt
#pipe2ucsf test.ft2 test.ucsf