#!/bin/csh

var2pipe -in ./fid \
 -noaswap  \
  -xN              3072  \
  -xT              1536  \
  -xMODE        Complex  \
  -xSW        17361.111  \
  -xOBS         799.910  \
  -xCAR           4.791  \
  -xLAB              H1  \
  -ndim               1  \
  -out ./test.fid -verb -ov

sleep 0
