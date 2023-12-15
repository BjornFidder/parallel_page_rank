#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --output=outlinks.out
##SBATCH -t 00:02:00

for n in 1 2 3 4 5 6 7 8 9 10
do
    for p in 1 4 8 16 32 64 128
    do
        # number of processor cores to use
        echo $p > input.tmp
        # vector length
        echo $N >> input.tmp

        srun "../bsp_pr" < "input.tmp"
        wait
    done
done