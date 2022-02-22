#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -ext -aswap -AMX -decim 1296 -dspfvs 20 -grpdly 67.9876556396484  \
  -xN              3072  -yN               460  \
  -xT              1484  -yT               230  \
  -xMODE            DQD  -yMODE    States-TPPI  \
  -xSW        15432.099  -ySW         5279.831  \
  -xOBS         700.303  -yOBS         176.116  \
  -xCAR           4.792  -yCAR         147.972  \
  -xLAB              1H  -yLAB             13C  \
  -ndim               2  -aq2D         Complex  \
  -out ./test.fid -verb -ov

sleep 0
