#!/bin/csh

bruk2pipe -in ./fid \
  -bad 0.0 -aswap -DMX -decim 1104 -dspfvs 20 -grpdly 67.9882354736328  \
  -xN              3840  \
  -xT              1810  \
  -xMODE            DQD  \
  -xSW        18115.942  \
  -xOBS         699.953  \
  -xCAR           4.791  \
  -xLAB              1H  \
  -ndim               1  \
  -out ./test.fid -verb -ov

sleep 5
