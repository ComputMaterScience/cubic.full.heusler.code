#!/bin/sh
#SBATCH -J test2 #job name
#SBATCH -p 52core #partition name
#SBATCH -N 1 #node
#SBATCH -n 52#core
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j

python ML_AHC_cal.py
