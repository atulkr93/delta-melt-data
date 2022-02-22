#!/bin/csh

bruk2pipe -in ./fid \
  -bad 0.0 -aswap -DMX -decim 1194.66666666667 -dspfvs 20 -grpdly 67.9878692626953  \
  -xN              3584  \
  -xT              1696  \
  -xMODE            DQD  \
  -xSW        16741.071  \
  -xOBS         700.303  \
  -xCAR           4.953  \
  -xLAB              1H  \
  -ndim               1  \
  -out ./test.fid -verb -ov

sleep 5
