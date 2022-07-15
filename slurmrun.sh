#!/bin/sh
#SBATCH --array=101-500
#SBATCH -c 1
#SBATCH --job-name=transport
#SBATCH --time=03:00:00
#SBATCH --mem=512
k=${SLURM_ARRAY_TASK_ID}
fname=runm800s33SI/positions${k}.dat
./cosmion.x 0.8 1.d-33 100000 ${fname} 0
