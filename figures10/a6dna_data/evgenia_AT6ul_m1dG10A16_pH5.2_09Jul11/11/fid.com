#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -noaswap -DMX -decim 24 -dspfvs 12 -grpdly 0  \
  -xN              1280  -yN               256  \
  -xT               573  -yT               128  \
  -xMODE            DQD  -yMODE  Echo-AntiEcho  \
  -xSW         7183.908  -ySW         3320.329  \
  -xOBS         600.133  -yOBS         150.925  \
  -xCAR           4.781  -yCAR         146.757  \
  -xLAB              1H  -yLAB             13C  \
  -ndim               2  -aq2D          States  \
  -out ./test.fid -verb -ov

sleep 5
