#!/bin/bash
 
# Request resources:
#SBATCH -c 1             # Number of CPU cores
#SBATCH --mem=32G         # memory required, <= 250G on standard nodes
#SBATCH --time=2-23:56:0 # time limit (format: days-hours:minutes:seconds)
#SBATCH --gres=tmp:12G   # temporary disk space ($TMPDIR), <= 400G
#SBATCH --job-name=%SCRIPTBASE%
#SBATCH --output=%SCRIPTBASE%.out
#SBATCH --error=%SCRIPTBASE%.err
# Run in the 'shared' queue (job may share node with other jobs)
#SBATCH -p shared
#SBATCH --export=ALL

# Commands to be run:

# Get current job's time limit in minutes, add 10 minute buffer
current_job=$SLURM_JOB_ID
time_limit=$(scontrol show job $current_job | grep TimeLimit | awk -F'TimeLimit=' '{print $2}' | cut -d' ' -f1)

# Parse time limit that could be either D-HH:MM:SS or HH:MM:SS
if [[ $time_limit == *-* ]]; then
    # Format with days: D-HH:MM:SS
    IFS=- read -r days rest <<< "$time_limit"
    IFS=: read -r hours minutes seconds <<< "$rest"
    delay_minutes=$(( days * 24 * 60 + hours * 60 + minutes + 10 ))
else
    # Format without days: HH:MM:SS
    IFS=: read -r hours minutes seconds <<< "$time_limit"
    delay_minutes=$(( hours * 60 + minutes + 10 ))
fi

if [[ $(squeue -u $USER --format "%j" | grep -c "__%SCRIPTBASE%") -gt 1 ]]; then
  echo "Another continuation job is queued. Exiting."
  exit 1
else
  echo "Queueing continuation job for " $delay_minutes " minutes."
  # Queue continuation job
  continuation_job=$(sbatch --parsable --begin=now+${delay_minutes}minutes --time 2:15:00 --job-name=__%SCRIPTBASE% --output=_%SCRIPTBASE%.out --error=_%SCRIPTBASE%.err slurm/%SCRIPTBASE%.sh)
fi

# Set git identifier
git config --global user.name "RevBayes serial analysis"

# Check out analysis repo
git clone --depth 1 https://${NT_WRITE}@github.com/neo-trans/%SCRIPTBASE%.git /nobackup/$USER/%SCRIPTBASE%
cd /nobackup/$USER/%SCRIPTBASE%
git checkout main
git pull origin main --rebase

# Run RevBayes
module load gcc/11.2 boost
start_time=$(date +%s)
../revbayes/projects/cmake/build/rb mc3serial.Rev --args 333

end_time=$(date +%s)
echo "RevBayes job terminated at $end_time."
elapsed=$((end_time - start_time))

if [ $elapsed -le 60 ]; then
  echo "RevBayes only ran for $elapsed seconds. Cancelling continuation job."
  scancel $continuation_job
fi

# Sync repo with any changes that have happened remotely during the run
git fetch --depth 1 origin

for file in *.trees; do
  tar -czf "${file%.trees}.tar.gz" "$file"
done

# Commit output files to git
git add %SCRIPTBASE%*.ckp  %SCRIPTBASE%*.log  %SCRIPTBASE%*.tar.gz
git commit -m "MCMCMC output files"
git rebase main
git push origin main
