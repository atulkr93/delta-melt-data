#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -aswap -DMX -decim 1296 -dspfvs 20 -grpdly 67.9876556396484  \
  -xN              2048  \
  -xT              1023  \
  -xMODE            DQD  \
  -xSW        15432.099  \
  -xOBS         700.178  \
  -xCAR           4.953  \
  -xLAB              1H  \
  -ndim               1  \
  -out ./test.fid -verb -ov

sleep 5
