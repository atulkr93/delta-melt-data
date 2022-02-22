#!/bin/csh

bruk2pipe -in ./ser \
  -bad 0.0 -ext -aswap -AMX -decim 5712 -dspfvs 20 -grpdly 67.9851226806641  \
  -xN               768  -yN               100  \
  -xT               349  -yT                50  \
  -xMODE            DQD  -yMODE  Echo-AntiEcho  \
  -xSW         3501.401  -ySW         1760.563  \
  -xOBS         700.104  -yOBS         176.055  \
  -xCAR           4.771  -yCAR          87.737  \
  -xLAB              1H  -yLAB             13C  \
  -ndim               2  -aq2D         Complex  \
  -out ./test.fid -verb -ov

sleep 0
