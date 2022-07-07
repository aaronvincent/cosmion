#!/bin/sh
#SBATCH -c 1
#SBATCH --time=24:00:00
#SBATCH --mem=2048
export OMP_NUM_THREADS=1
k=${SLURM_ARRAY_TASK_ID}
fname=positions${k}.dat
./cosmion.x 5. 1.d-37 10000 fname
