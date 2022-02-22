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
| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 1 -c 0.5    \
| nmrPipe  -fn ZF -size 8192                               \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 214.70 -p1 -71.00 -di -verb         \
| nmrPipe  -fn EXT -xn 8.4ppm -x1 6.8ppm -sw -verb    \
| nmrPipe  -fn TP                                     \
| nmrPipe  -fn SP -off 0.5 -end 1.00 -pow 1 -c 0.5    \
| nmrPipe  -fn ZF -size 4096                               \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 101.70 -p1 -46.00 -di -verb         \
| nmrPipe  -fn EXT -xn 156ppm -x1 135ppm -sw -verb    \
   -ov -out test.ft2
#pipe2ucsf test.ft2 scaf2-TGC_low_5.4_HSQC_30C_Aro.ucsf
