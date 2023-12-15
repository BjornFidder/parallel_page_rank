#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --output=outlinks.out
##SBATCH -t 00:02:00

for n in 1 2 3 4 5 6 7 8 9 10
do
    N=$((n*1000000))
    echo $N > input.tmp
    srun "../seq_pr" < "input.tmp"
    wait

    for p in 1 4 8 16 32 64 128
    do
        # BSPedupack programs read arguments from stdin, so we have to "pipe them in"
        # number of processor cores to use
        echo $p > input.tmp
        # vector length
        echo $N >> input.tmp

        srun "../bsp_pr" < "input.tmp"
        wait
    done
done