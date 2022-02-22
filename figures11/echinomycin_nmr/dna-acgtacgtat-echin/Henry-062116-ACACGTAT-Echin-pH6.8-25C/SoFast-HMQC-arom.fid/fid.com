#!/bin/csh

var2pipe -in ./fid \
 -noaswap  \
  -xN              2048  -yN               360  \
  -xT              1024  -yT               180  \
  -xMODE        Complex  -yMODE        Complex  \
  -xSW        16025.641  -ySW         4425.637  \
  -xOBS         799.910  -yOBS         201.165  \
  -xCAR           4.791  -yCAR         145.150  \
  -xLAB              H1  -yLAB             C13  \
  -ndim               2  -aq2D          States  \
  -out ./test.fid -verb -ov

sleep 0
