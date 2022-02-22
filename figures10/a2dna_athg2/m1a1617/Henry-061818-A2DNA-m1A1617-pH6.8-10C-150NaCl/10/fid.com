#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -aswap -DMX -decim 2600 -dspfvs 20 -grpdly 67.9858245849609  \
  -xN              1536  -yN               150  \
  -xT               761  -yT                75  \
  -xMODE            DQD  -yMODE    States-TPPI  \
  -xSW         7692.308  -ySW         2465.483  \
  -xOBS         700.178  -yOBS         176.084  \
  -xCAR           4.953  -yCAR         142.642  \
  -xLAB              HN  -yLAB             13C  \
  -ndim               2  -aq2D          States  \
  -out ./test.fid -verb -ov

sleep 5
