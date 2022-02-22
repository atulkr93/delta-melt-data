#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -aswap -DMX -decim 1792 -dspfvs 20 -grpdly 67.9841766357422  \
  -xN              2048  -yN               440  \
  -xT              1024  -yT               220  \
  -xMODE            DQD  -yMODE  Echo-AntiEcho  \
  -xSW        11160.714  -ySW         4930.966  \
  -xOBS         700.303  -yOBS         176.116  \
  -xCAR           4.791  -yCAR         144.750  \
  -xLAB              1H  -yLAB             13C  \
  -ndim               2  -aq2D          States  \
  -out ./test.fid -verb -ov

sleep 5
