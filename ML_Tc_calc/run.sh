#!/bin/sh
#SBATCH -J test #job name
#SBATCH -p 52core #partition name
#SBATCH -N 1 #node
#SBATCH -n 52 #core
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j

export OMP_NUM_THREADS=1 #number of threads

python ML_Tc_cal.py
