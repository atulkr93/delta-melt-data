#!/bin/csh

nmrPipe -in test.fid \
#| nmrPipe  -fn POLY -time \
#| nmrPipe  -fn SOL                           \
| nmrPipe  -fn GM -g1 5 -g2 15 -g3 0.0 -c 1.0 \
#| nmrPipe  -fn SP -off 0.45 -end 0.98 -pow 1 -c 0.5   \
| nmrPipe  -fn ZF -size 8192                               \
| nmrPipe  -fn FT                                     \
| nmrPipe  -fn PS -p0 -90  -p1 -180 -di                 \
| nmrPipe  -fn EXT -x1 15.0ppm -xn 12.0ppm -sw    \
| nmrPipe  -fn TP                                     \
#| nmrPipe  -fn GM -g1 5 -g2 15 -g3 0.0 -c 1.0 \
| nmrPipe  -fn SP -off 0.45 -end 0.98 -pow 1 -c 0.5   \
| nmrPipe  -fn ZF -size 2048                               \
| nmrPipe  -fn FT -auto                               \
#| nmrPipe -fn REV -sw      \
| nmrPipe  -fn PS -p0 -90 -p1 180.0 -di                \
#| nmrPipe  -fn EXT -x1 139.0ppm -xn 144.0ppm -sw \
   -verb -ov -out test.ft2
pipe2ucsf test.ft2 A6DNAm3T9_sofast_HMQC_Imino.ucsf
