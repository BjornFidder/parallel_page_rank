#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
##SBATCH -t 00:10:00

P=1

# BSPedupack programs read arguments from stdin, so we have to "pipe them in"
# number of processor cores to use
echo $P > input.tmp

srun "../bench" < "input.tmp"
wait
