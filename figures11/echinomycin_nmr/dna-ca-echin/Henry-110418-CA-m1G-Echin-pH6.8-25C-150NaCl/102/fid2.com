#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -ext -aswap -AMX -decim 1296 -dspfvs 20 -grpdly 67.9876556396484  \
  -xN              3072  -yN               144  \
  -xT              1484  -yT                72  \
  -xMODE            DQD  -yMODE    States-TPPI  \
  -xSW        15432.099  -ySW         1760.563  \
  -xOBS         700.178  -yOBS         176.084  \
  -xCAR           4.791  -yCAR         141.691  \
  -xLAB              1H  -yLAB             13C  \
  -ndim               2  -aq2D         Complex  \
  -out ./test.fid -verb -ov

sleep 0
