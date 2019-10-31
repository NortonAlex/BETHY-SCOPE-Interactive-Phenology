#!/bin/bash
#PBS -P w22
#PBS -q copyq
#PBS -l walltime=00:05:00
#PBS -l mem=100MB 
#PBS -l ncpus=1
#PBS -l wd

python write_block_control_files.py
