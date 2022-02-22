#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -ext -aswap -AMX -decim 1296 -dspfvs 20 -grpdly 67.9876556396484  \
  -xN              4096  \
  -xT              2048  \
  -xMODE            DQD  \
  -xSW        15432.099  \
  -xOBS         700.066  \
  -xCAR           4.755  \
  -xLAB              HN  \
  -ndim               1  \
  -out ./test.fid -verb -ov

sleep 0
