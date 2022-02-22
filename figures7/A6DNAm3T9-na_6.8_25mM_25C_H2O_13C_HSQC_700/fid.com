#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -aswap -DMX -decim 1792 -dspfvs 20 -grpdly 67.9841766357422  \
  -xN              2048  -yN               256  \
  -xT              1024  -yT               128  \
  -xMODE            DQD  -yMODE  Echo-AntiEcho  \
  -xSW        11160.714  -ySW         4930.966  \
  -xOBS         700.178  -yOBS         176.084  \
  -xCAR           4.771  -yCAR         144.735  \
  -xLAB              1H  -yLAB             13C  \
  -ndim               2  -aq2D          States  \
  -out ./test.fid -verb -ov

sleep 0
