#/bin/bash
#PBS -P w22
#PBS -q normal
#PBS -l walltime=00:40:00
#PBS -l mem=32GB 
#PBS -l ncpus=16
#PBS -l wd

fluoro_folder=/home/563/ajn563/short/fluoro2
cd $fluoro_folder

export OMP_NUM_THREADS=16
## lift any arbitrary limits on the per-process stack-size
ulimit -s unlimited
## raise the limit for the OpenMP per-thread stack-size
export OMP_STACKSIZE=100MB

#module load netcdf
#!make post

./post

