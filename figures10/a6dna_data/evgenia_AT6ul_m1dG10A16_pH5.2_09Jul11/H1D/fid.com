#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -ext -noaswap -AMX -decim 12 -dspfvs 12 -grpdly 0  \
  -xN              3840  \
  -xT              1800  \
  -xMODE            DQD  \
  -xSW        15015.015  \
  -xOBS         600.133  \
  -xCAR           4.781  \
  -xLAB              1H  \
  -ndim               1  \
  -out ./test.fid -verb -ov

sleep 0
