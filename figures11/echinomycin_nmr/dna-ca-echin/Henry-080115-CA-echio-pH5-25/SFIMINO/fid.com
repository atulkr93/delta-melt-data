#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -aswap -DMX -decim 1144 -dspfvs 20 -grpdly 67.9914703369141  \
  -xN              4352  \
  -xT              2096  \
  -xMODE            DQD  \
  -xSW        17482.517  \
  -xOBS         699.953  \
  -xCAR           4.791  \
  -xLAB              1H  \
  -ndim               1  \
  -out ./test.fid -verb -ov

sleep 5
