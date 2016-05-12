#/bin/bash
#PBS -P w22
#PBS -q express
#PBS -l walltime=01:00:00
#PBS -l mem=8GB 
#PBS -l ncpus=16
#PBS -l wd

fluoro_folder=/home/563/ajn563/short/fluoro2
cd $fluoro_folder

export OMP_NUM_THREADS=16
## lift any arbitrary limits on the per-process stack-size
ulimit -s unlimited
## raise the limit for the OpenMP per-thread stack-size
export OMP_STACKSIZE=100MB

module load netcdf
#!make post

./post

