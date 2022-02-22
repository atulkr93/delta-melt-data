#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -aswap -DMX -decim 1904 -dspfvs 20 -grpdly 67.9868774414062  \
  -xN              2048  -yN               290  \
  -xT              1024  -yT               145  \
  -xMODE            DQD  -yMODE  Echo-AntiEcho  \
  -xSW        10504.202  -ySW         5813.953  \
  -xOBS         700.178  -yOBS         176.072  \
  -xCAR           4.771  -yCAR          77.235  \
  -xLAB              1H  -yLAB             13C  \
  -ndim               2  -aq2D          States  \
  -out ./test.fid -verb -ov

sleep 0
