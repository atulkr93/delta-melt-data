#!/bin/csh

#
# Basic 2D Phase-Sensitive Processing:
#   Cosine-Bells are used in both dimensions.
#   Use of "ZF -auto" doubles size, then rounds to power of 2.
#   Use of "FT -auto" chooses correct Transform mode.
#   Imaginaries are deleted with "-di" in each dimension.
#   Phase corrections should be inserted by hand.

nmrPipe -in test.fid \
| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 1 -c 0.5    \
| nmrPipe  -fn ZF -size 8192                               \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 11.20 -p1 216.00 -di -verb         \
| nmrPipe  -fn EXT -x1 6.6ppm -xn 8.4ppm -sw \
| nmrPipe  -fn TP                                     \
| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 1 -c 0.5    \
| nmrPipe  -fn ZF -size 2048                               \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 -58.00 -p1 -72.00 -di -verb         \
| nmrPipe  -fn EXT -x1 146ppm -xn 137ppm -sw \
   -ov -out test.ft2
#pipe2ucsf test.ft2 S1lb_rA16dsA6DNA_pH5.4_25C_Caro.ucsf
