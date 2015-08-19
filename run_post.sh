#!/bin/bash

export OMP_NUM_THREADS=16
## lift any arbitrary limits on the per-process stack-size
ulimit -s unlimited
## raise the limit for the OpenMP per-thread stack-size
export OMP_STACKSIZE=100MB

./post

