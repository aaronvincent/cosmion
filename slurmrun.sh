#!/bin/sh
#SBATCH --array=0-20
#SBATCH -c 1
#SBATCH --job-name=transport
#SBATCH --time=03:00:00
#SBATCH --mem=1024
k=${SLURM_ARRAY_TASK_ID}
fname=run1/positions${k}.dat
./cosmion.x 5. 1.d-37 100000 ${fname}
