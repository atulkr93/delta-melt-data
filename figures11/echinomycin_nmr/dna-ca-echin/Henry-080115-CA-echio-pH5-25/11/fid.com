#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -aswap -DMX -decim 1144 -dspfvs 20 -grpdly 67.9914703369141  \
  -xN              4096  -yN              1600  \
  -xT              2048  -yT               800  \
  -xMODE            DQD  -yMODE    States-TPPI  \
  -xSW        17482.517  -ySW        17482.517  \
  -xOBS         699.953  -yOBS         699.953  \
  -xCAR           4.791  -yCAR           4.791  \
  -xLAB             1Hx  -yLAB             1Hy  \
  -ndim               2  -aq2D          States  \
  -out ./test.fid -verb -ov

sleep 5
