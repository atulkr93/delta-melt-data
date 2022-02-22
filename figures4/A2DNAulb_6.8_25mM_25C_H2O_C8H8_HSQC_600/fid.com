#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -ext -aswap -AMX -decim 3328 -dspfvs 20 -grpdly 67.9842681884766  \
  -xN              1024  -yN               200  \
  -xT               512  -yT               100  \
  -xMODE            DQD  -yMODE  Echo-AntiEcho  \
  -xSW         6009.615  -ySW         3770.739  \
  -xOBS         599.663  -yOBS         150.806  \
  -xCAR           4.771  -yCAR         144.719  \
  -xLAB              1H  -yLAB             13C  \
  -ndim               2  -aq2D         Complex  \
  -out ./test.fid -verb -ov

sleep 0
