#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -ext -aswap -AMX -decim 4170.66666666667 -dspfvs 20 -grpdly 67.9881591796875  \
  -xN              1024  -yN               120  \
  -xT               479  -yT                60  \
  -xMODE            DQD  -yMODE  Echo-AntiEcho  \
  -xSW         4795.396  -ySW         1760.563  \
  -xOBS         599.663  -yOBS         150.797  \
  -xCAR           4.771  -yCAR          87.733  \
  -xLAB              1H  -yLAB             13C  \
  -ndim               2  -aq2D         Complex  \
  -out ./test.fid -verb -ov

sleep 0
