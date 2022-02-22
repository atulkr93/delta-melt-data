#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -aswap -DMX -decim 1392 -dspfvs 20 -grpdly 67.9874420166016  \
  -xN              3328  -yN              2264  \
  -xT              1561  -yT              1132  \
  -xMODE            DQD  -yMODE    States-TPPI  \
  -xSW        14367.816  -ySW        14409.222  \
  -xOBS         599.663  -yOBS         599.663  \
  -xCAR           4.791  -yCAR           4.791  \
  -xLAB             1Hx  -yLAB             1Hy  \
  -ndim               2  -aq2D          States  \
  -out ./test.fid -verb -ov

sleep 5
