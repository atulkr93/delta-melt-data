#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -aswap -DMX -decim 1520 -dspfvs 20 -grpdly 67.9868469238281  \
  -xN              3072  -yN              1280  \
  -xT              1536  -yT               640  \
  -xMODE            DQD  -yMODE    States-TPPI  \
  -xSW        13157.895  -ySW        12285.012  \
  -xOBS         599.663  -yOBS         599.663  \
  -xCAR           4.771  -yCAR           4.771  \
  -xLAB             1Hx  -yLAB             1Hy  \
  -ndim               2  -aq2D          States  \
  -out ./test.fid -verb -ov

sleep 0
