#!/bin/csh

nmrPipe -in test.fid \
#| nmrPipe  -fn POLY -time \
#| nmrPipe  -fn SOL                           \
| nmrPipe  -fn GM -g1 5 -g2 10 -g3 0.0	-c 1\
#| nmrPipe  -fn SP -off 0.45 -end 0.98 -pow 1 -c 0.5    \
| nmrPipe  -fn ZF -size 4096                               \
| nmrPipe  -fn FT                                     \
| nmrPipe  -fn PS -p0 70.0  -p1 0.0 -di                 \
#| nmrPipe  -fn POLY -auto \
| nmrPipe  -fn EXT -x1 6.5ppm -xn 5.0ppm -sw    \
| nmrPipe  -fn TP                                     \
#| nmrPipe  -fn GM -g1 5 -g2 15 -g3 0.0 -c 1.0 \
| nmrPipe  -fn SP -off 0.55 -end 0.98 -pow 2 -c 0.5    \
| nmrPipe  -fn ZF -size 2048                               \
| nmrPipe  -fn FT -auto                               \
| nmrPipe  -fn PS -p0 90 -p1 0.0 -di                \
| nmrPipe  -fn EXT -x1 90.0ppm -xn 80.0ppm -sw \
   -verb -ov -out test.ft2
#pipe2ucsf test.ft2  A6DNAm1A16_T6lb_HSQC_C1p.ucsf
