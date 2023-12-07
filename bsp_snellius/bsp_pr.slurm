#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
##SBATCH -t 00:01:00

N=100
echo $N > input.tmp
srun "../seq_pr" < "input.tmp"
wait

for p in 1 4 8 16
do
    # BSPedupack programs read arguments from stdin, so we have to "pipe them in"
    # number of processor cores to use
    echo $p > input.tmp
    # vector length
    echo $N >> input.tmp

    # runs the inner product (inprod) example using 4 processes
    srun "../bsp_pr" < "input.tmp"
    wait
done