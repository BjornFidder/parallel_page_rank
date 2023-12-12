#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
##SBATCH -t 00:01:00

for p in 1 4 8 16 32 64 128
do
    # BSPedupack programs read arguments from stdin, so we have to "pipe them in"
    # number of processor cores to use
    echo $p > input.tmp

    # runs the inner product (inprod) example using 4 processes
    srun "../bench" < "input.tmp"
    wait
done