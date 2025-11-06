#!/bin/bash
 
# Request resources:
#SBATCH -n 16          # number of tasks, up to 128. ## Replacing "SBATCH -c 16" : # number of CPU cores, one per thread, up to 128
#SBATCH --mem=24G      # memory required, up to 250G on standard nodes
#SBATCH --time=2-23:45:0  # time limit for job (format:  days-hours:minutes:seconds)
#SBATCH --gres=tmp:128G # temporary disk space required on the compute node ($TMPDIR), up to 400G
#SBATCH --job-name=%SCRIPTBASE%
#SBATCH --output=/nobackup/%u/neotrans/%SCRIPTBASE%.out
#SBATCH --error=/nobackup/%u/neotrans/%SCRIPTBASE%.err
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
module load gcc/11.2 boost openmpi/4.1.1
start_time=$(date +%s)
mpirun /nobackup/$USER/revbayes/projects/cmake/build-mpi/rb-mpi mcmcmc.Rev --args 333

end_time=$(date +%s)
echo "RevBayes job terminated at $end_time."

# Restore original .Rev script to avoid conflicts
git restore mcmcmc.Rev

# Sync repo with any changes that have happened remotely during the run
git fetch --depth 1 origin

if ls *.trees 1> /dev/null 2>&1; then
  for file in *.trees; do
    tar -czf "${file%.trees}.tar.gz" "$file"
  done
  
  # Commit output files to git
  git add %SCRIPTBASE%*.ckp %SCRIPTBASE%*.log %SCRIPTBASE%*.var %SCRIPTBASE%*.tar.gz
  git commit -m "MCMCMC output files"
  git rebase main
  git push origin main
else 
  echo "No .trees files found. No output to commit."
fi

# Record usage of tmpDir
du -hs $TMPDIR > /nobackup/$USER/%SCRIPTBASE%/mc3-tmpdir_usage.log
