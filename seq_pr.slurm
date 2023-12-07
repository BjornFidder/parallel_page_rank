#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
##SBATCH -t 00:01:00

# BSPedupack programs read arguments from stdin, so we have to "pipe them in"
# vector length
echo 100000 >> input.tmp

# runs the sequential pagerank
srun seq_pr < "input.tmp"