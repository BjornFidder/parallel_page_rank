#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --output=put.out
#SBATCH -t 00:10:00

for P in 1 2 4 8 16 32 64 128 
do
    # BSPedupack programs read arguments from stdin, so we have to "pipe them in"
    # number of processor cores to use
    echo $P > input.tmp

    srun "../bench" < "input.tmp"
    wait
done


