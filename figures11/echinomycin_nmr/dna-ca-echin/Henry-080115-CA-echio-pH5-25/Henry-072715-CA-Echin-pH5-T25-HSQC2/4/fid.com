#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -aswap -DMX -decim 1520 -dspfvs 20 -grpdly 67.9868469238281  \
  -xN              2816  -yN               472  \
  -xT              1300  -yT               236  \
  -xMODE            DQD  -yMODE    States-TPPI  \
  -xSW        13157.895  -ySW         3921.569  \
  -xOBS         599.663  -yOBS         150.806  \
  -xCAR           4.791  -yCAR         145.261  \
  -xLAB              1H  -yLAB             13C  \
  -ndim               2  -aq2D          States  \
  -out ./test.fid -verb -ov

sleep 0
