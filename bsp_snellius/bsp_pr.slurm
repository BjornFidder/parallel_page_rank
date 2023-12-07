#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
##SBATCH -t 00:01:00

# BSPedupack programs read arguments from stdin, so we have to "pipe them in"
# number of processor cores to use
echo 4 > input.tmp
# vector length
echo 100000 >> input.tmp

# runs the inner product (inprod) example using 4 processes
srun "../bsp_pr" < "input.tmp"