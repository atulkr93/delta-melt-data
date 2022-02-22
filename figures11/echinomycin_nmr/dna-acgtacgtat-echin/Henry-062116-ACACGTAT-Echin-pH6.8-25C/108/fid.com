#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -aswap -DMX -decim 1296 -dspfvs 20 -grpdly 67.9876556396484  \
  -xN              3072  -yN               110  \
  -xT              1536  -yT                55  \
  -xMODE            DQD  -yMODE    States-TPPI  \
  -xSW        15432.099  -ySW         1561.524  \
  -xOBS         700.303  -yOBS          70.972  \
  -xCAR           4.791  -yCAR         155.089  \
  -xLAB              HN  -yLAB             15N  \
  -ndim               2  -aq2D          States  \
  -out ./test.fid -verb -ov

sleep 5
