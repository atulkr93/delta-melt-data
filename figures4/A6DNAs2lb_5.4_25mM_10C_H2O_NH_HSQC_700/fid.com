#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -aswap -DMX -decim 1296 -dspfvs 20 -grpdly 67.9876556396484  \
  -xN              3328  -yN               150  \
  -xT              1542  -yT                75  \
  -xMODE            DQD  -yMODE    States-TPPI  \
  -xSW        15432.099  -ySW         1844.338  \
  -xOBS         700.178  -yOBS          70.959  \
  -xCAR           4.915  -yCAR         147.773  \
  -xLAB              HN  -yLAB             15N  \
  -ndim               2  -aq2D          States  \
  -out ./test.fid -verb -ov

sleep 0
