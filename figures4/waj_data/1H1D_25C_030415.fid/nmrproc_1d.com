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
| nmrPipe  -fn ZF -size 16384                          \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn EXT -x1 14.3ppm -xn 12.5ppm -sw \
| nmrPipe  -fn PS -p0 18.00 -p1 0.00 -di -verb         \
| nmrPipe  -fn POLY -auto \
   -ov -out test.ft2
