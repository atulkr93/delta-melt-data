#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -ext -aswap -AMX -decim 2088 -dspfvs 20 -grpdly 67.9876556396484  \
  -xN              2048  -yN                 6  \
  -xT               957  -yT                 6  \
  -xMODE            DQD  -yMODE           Real  \
  -xSW         9578.544  -ySW            6.000  \
  -xOBS         599.663  -yOBS           1.000  \
  -xCAR           4.867  -yCAR           0.000  \
  -xLAB              1H  -yLAB             TAU  \
  -ndim               2  -aq2D         Complex  \
  -out ./test.fid -verb -ov

sleep 0
