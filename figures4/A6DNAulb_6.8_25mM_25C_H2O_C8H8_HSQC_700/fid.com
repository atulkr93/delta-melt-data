#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -aswap -DMX -decim 2784 -dspfvs 20 -grpdly 67.9842529296875  \
  -xN              1536  -yN               250  \
  -xT               718  -yT               125  \
  -xMODE            DQD  -yMODE  Echo-AntiEcho  \
  -xSW         7183.908  -ySW         4401.408  \
  -xOBS         700.178  -yOBS         176.084  \
  -xCAR           4.771  -yCAR         145.236  \
  -xLAB              1H  -yLAB             13C  \
  -ndim               2  -aq2D          States  \
  -out ./test.fid -verb -ov

sleep 0
