#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -ext -aswap -AMX -decim 2784 -dspfvs 20 -grpdly 67.9842529296875  \
  -xN              1536  -yN               512  \
  -xT               717  -yT               256  \
  -xMODE            DQD  -yMODE  Echo-AntiEcho  \
  -xSW         7183.908  -ySW         4401.408  \
  -xOBS         699.953  -yOBS         176.028  \
  -xCAR           4.790  -yCAR         145.165  \
  -xLAB              1H  -yLAB             13C  \
  -ndim               2  -aq2D         Complex  \
  -out ./test.fid -verb -ov

sleep 0
