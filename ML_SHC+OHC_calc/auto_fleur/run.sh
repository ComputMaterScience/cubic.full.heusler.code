#!/bin/sh
#SBATCH -J test #job name
#SBATCH -p 36core #partition name
#SBATCH -N 1 #node
#SBATCH -n 36 #core
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j

export OMP_NUM_THREADS=36 #number of threads

python auto_shc_wannier90_fleur.py
