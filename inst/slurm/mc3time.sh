#!/bin/bash
 
# Request resources:
#SBATCH -n 1           # number of tasks, up to 128.
#SBATCH --mem=24G      # memory required, up to 250G on standard nodes
#SBATCH --time=0-12:34:56  # time limit for job (format:  days-hours:minutes:seconds)
#SBATCH --gres=tmp:128G # temporary disk space required on the compute node ($TMPDIR), up to 400G
#SBATCH --job-name=time-%SCRIPTBASE%
#SBATCH --output=/nobackup/%u/neotrans/time-%SCRIPTBASE%.out
#SBATCH --error=/nobackup/%u/neotrans/time-%SCRIPTBASE%.err
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

# Set RevBayes run time
sed -i 's/runHours = [0-9]\+/runHours = 71/' mcmcmc.Rev

# Run RevBayes
module load gcc/11.2 boost
start_time=$(date +%s)

/usr/bin/time -v -o /nobackup/$USER/%SCRIPTBASE%/%SCRIPTBASE%.mem.log \
  /nobackup/$USER/revbayes/projects/cmake/build/rb mc3time.Rev

end_time=$(date +%s)
echo "RevBayes job terminated at $end_time."

# Record usage of tmpDir
du -hs $TMPDIR > /nobackup/$USER/%SCRIPTBASE%/%SCRIPTBASE%.tmp_usage.log

# Restore original .Rev script to avoid conflicts
git restore mc3time.Rev

# Sync repo with any changes that have happened remotely during the run
git fetch --depth 1 origin

if ls *.trees 1> /dev/null 2>&1; then
  for file in *.trees; do
    tar -czf "${file%.trees}.tar.gz" "$file"
  done
  
  # Commit output files to git
  git add %SCRIPTBASE%*.ckp  %SCRIPTBASE%*.log  %SCRIPTBASE%*.tar.gz mc3time.Rev
  git commit -m "MCMCMC timing files"
  git pull --rebase
  git push origin main
else 
  echo "No .trees files found. No output to commit."
fi

