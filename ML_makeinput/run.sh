#!/bin/sh
#SBATCH -J test #job name
#SBATCH -p node #partition name
#SBATCH -N 1 #node
#SBATCH -n 4 #core
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so #Do not change here!!

/home/user2/thiho/MLdata/getdata
