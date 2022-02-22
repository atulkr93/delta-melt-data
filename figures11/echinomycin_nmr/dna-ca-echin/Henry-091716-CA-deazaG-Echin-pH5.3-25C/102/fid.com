#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -aswap -DMX -decim 1296 -dspfvs 20 -grpdly 67.9876556396484  \
  -xN              3072  -yN               460  \
  -xT              1484  -yT               230  \
  -xMODE            DQD  -yMODE    States-TPPI  \
  -xSW        15432.099  -ySW         5279.831  \
  -xOBS         700.303  -yOBS         176.090  \
  -xCAR           4.853  -yCAR         147.633  \
  -xLAB              1H  -yLAB             13C  \
  -ndim               2  -aq2D          States  \
  -out ./test.fid -verb -ov

sleep 0
