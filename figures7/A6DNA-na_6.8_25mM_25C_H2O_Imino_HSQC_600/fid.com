#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -aswap -DMX -decim 1520 -dspfvs 20 -grpdly 67.9868469238281  \
  -xN              2816  -yN                88  \
  -xT              1315  -yT                44  \
  -xMODE            DQD  -yMODE    States-TPPI  \
  -xSW        13157.895  -ySW         1215.362  \
  -xOBS         599.663  -yOBS          60.772  \
  -xCAR           4.773  -yCAR         152.073  \
  -xLAB              HN  -yLAB             15N  \
  -ndim               2  -aq2D          States  \
  -out ./test.fid -verb -ov

sleep 0
