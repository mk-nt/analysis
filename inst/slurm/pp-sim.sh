#!/bin/bash
 
# Request resources:
#SBATCH -n 1           # number of tasks, up to 128.
#SBATCH --mem=24G      # memory required, up to 250G on standard nodes
#SBATCH --time=0-21:34:56  # time limit for job (format:  days-hours:minutes:seconds)
#SBATCH --gres=tmp:128G # temporary disk space required on the compute node ($TMPDIR), up to 400G
#SBATCH --job-name=%STAT%-%SCRIPTBASE%
#SBATCH --output=/nobackup/%u/%SCRIPTBASE%/%STAT%-%SCRIPTBASE%.out
#SBATCH --error=/nobackup/%u/neotrans/%STAT%-%SCRIPTBASE%.err
#SBATCH -p shared
#SBATCH --export=ALL

# Set git identifier
git config --global user.name "RevBayes analysis"

# Check out analysis repo - clone if doesn't exist, pull if it does
if [ -d "/nobackup/$USER/%SCRIPTBASE%/.git" ]; then
  echo "Repository already exists, pulling latest changes..."
  cd /nobackup/$USER/%SCRIPTBASE%
  git checkout main
  git pull origin main --rebase
else
  echo "Cloning repository..."
  git clone --depth 1 https://${NT_WRITE}@github.com/neo-trans/%SCRIPTBASE%.git /nobackup/$USER/%SCRIPTBASE%
  cd /nobackup/$USER/%SCRIPTBASE%
  git checkout main
fi

# Run RevBayes
module load gcc/11.2 boost
start_time=$(date +%s)

/nobackup/$USER/revbayes/projects/cmake/build/rb \
  ppsim_%SCRIPTID%.Rev %PID% %BURNIN% %STAT%

end_time=$(date +%s)
echo "RevBayes job terminated at $end_time."

# Sync repo with any changes that have happened remotely during the run
git fetch --depth 1 origin

# Commit output files to git
git add %STAT%-%SCRIPTBASE%.out
git commit -m "Simulate datasets from posterior"
git pull --rebase
git push origin main
