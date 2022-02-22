#!/bin/csh

bruk2pipe -in ./fid \
  -bad 0.0 -aswap -DMX -decim 1296 -dspfvs 20 -grpdly 67.9876556396484  \
  -xN              4864  \
  -xT              2314  \
  -xMODE            DQD  \
  -xSW        15432.099  \
  -xOBS         599.663  \
  -xCAR           4.791  \
  -xLAB              1H  \
  -ndim               1  \
  -out ./test.fid -verb -ov

sleep 5
