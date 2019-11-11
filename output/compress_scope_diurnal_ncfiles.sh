#!/bin/bash
#PBS -P w22
#PBS -q copyq
#PBS -l walltime=00:10:00
#PBS -l mem=1GB 
#PBS -l ncpus=1
#PBS -l wd

## Name the parameter experiment you want to compress 
#parname=1000

python compress_scope_diurnal_ncfiles.py  #${parname}
