#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -ext -aswap -AMX -decim 2784 -dspfvs 20 -grpdly 67.9842529296875  \
  -xN              1536  -yN               200  \
  -xT               717  -yT               100  \
  -xMODE            DQD  -yMODE  Echo-AntiEcho  \
  -xSW         7183.908  -ySW         2111.486  \
  -xOBS         599.663  -yOBS         150.806  \
  -xCAR           4.953  -yCAR         144.925  \
  -xLAB              1H  -yLAB             13C  \
  -ndim               2  -aq2D         Complex  \
  -out ./test.fid -verb -ov

sleep 0
