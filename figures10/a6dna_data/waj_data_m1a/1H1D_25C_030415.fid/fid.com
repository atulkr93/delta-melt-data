#!/bin/csh

var2pipe -in ./fid \
 -noaswap  \
  -xN              3352  \
  -xT              1676  \
  -xMODE        Complex  \
  -xSW        18939.394  \
  -xOBS         799.910  \
  -xCAR           4.791  \
  -xLAB              H1  \
  -ndim               1  \
  -out ./test.fid -verb -ov

sleep 0
