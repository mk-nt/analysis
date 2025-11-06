#!/bin/bash
 
#SBATCH --array=1-4
#SBATCH --job-name=ts-%PID%
#SBATCH --output=/nobackup/%u/neotrans/ts-%PID%_%a.out
#SBATCH --error=/nobackup/%u/neotrans/ts-%PID%_%a.err

#SBATCH -n 1            # number of tasks, up to 128.
#SBATCH --time=2-23:0:0    # time limit for job (format:  days-hours:minutes:seconds)
#SBATCH --mem=1600M     # memory required, up to 250G on standard nodes
#SBATCH --gres=tmp:4G   # temporary disk space required on $TMPDIR, up to 400G
#SBATCH -p shared

git config user.name "TreeSearch analysis"

cd /nobackup/$USER/%SCRIPTBASE%
git checkout main
git pull origin main --rebase

k_values=("1" "3" "10" "Inf")
array_index=$SLURM_ARRAY_TASK_ID
k_value=${k_values[$((array_index - 1))]}

module load r/4.1.2
Rscript "../neotrans/rScripts/TreeSearch.R" "$k_value"

# Record usage of tmpDir
du -hs $TMPDIR > /nobackup/$USER/%SCRIPTBASE%/ts-tmpdir_usage.log
