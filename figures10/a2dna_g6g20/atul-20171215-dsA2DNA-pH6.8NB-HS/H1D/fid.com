#!/bin/csh

bruk2pipe -in ./fid \
  -bad 0.0 -aswap -DMX -decim 1296 -dspfvs 20 -grpdly 67.9876556396484  \
  -xN              3328  \
  -xT              1542  \
  -xMODE            DQD  \
  -xSW        15432.099  \
  -xOBS         700.178  \
  -xCAR           4.791  \
  -xLAB              1H  \
  -ndim               1  \
  -out ./test.fid -verb -ov

sleep 0
