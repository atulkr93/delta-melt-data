#!/bin/csh

#
# Basic 2D Phase-Sensitive Processing:
#   Cosine-Bells are used in both dimensions.
#   Use of "ZF -auto" doubles size, then rounds to power of 2.
#   Use of "FT -auto" chooses correct Transform mode.
#   Imaginaries are deleted with "-di" in each dimension.
#   Phase corrections should be inserted by hand.

nmrPipe -in test.fid \
#| nmrPipe  -fn GM -g1 3 -g2 20 -g3 0.0 -c 1.0 \
#| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 1 -c 0.5    \
| nmrPipe  -fn ZF -size 16384                               \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 72.00 -p1 0.00 -di -verb         \
| nmrPipe  -fn EXT -x1 11ppm -xn 15ppm -sw \
| nmrPipe  -fn TP                                     \
#| nmrPipe  -fn GM -g1 5 -g2 15 -g3 0.0 -c 1.0 \
| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 1 -c 0.5    \
| nmrPipe  -fn ZF -size 4096                               \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 90 -p1 180 -di -verb         \
   -ov -out test.ft2
#pipe2ucsf test.ft2 A6DNA_sofast_HMQC_Imino.ucsf

