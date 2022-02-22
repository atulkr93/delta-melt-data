#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -ext -aswap -AMX -decim 1520 -dspfvs 20 -grpdly 67.9868469238281  \
  -xN              2048  \
  -xT              1023  \
  -xMODE            DQD  \
  -xSW        13157.895  \
  -xOBS         599.663  \
  -xCAR           4.791  \
  -xLAB              HN  \
  -ndim               1  \
  -out ./test.fid -verb -ov

sleep 0
