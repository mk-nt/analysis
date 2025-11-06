#!/bin/bash
 
# Request resources:
#SBATCH -n 16          # number of tasks, up to 128. ## Replacing "SBATCH -c 16" : # number of CPU cores, one per thread, up to 128
#SBATCH --mem=64G      # memory required, up to 250G on standard nodes
#SBATCH --time=2-23:56:0  # time limit for job (format:  days-hours:minutes:seconds)
#SBATCH --gres=tmp:128G # temporary disk space required on the compute node ($TMPDIR), up to 400G
#SBATCH --job-name=ml-%SCRIPTBASE%
#SBATCH --output=/nobackup/%u/neotrans/ml-%SCRIPTBASE%.out
#SBATCH --error=/nobackup/%u/neotrans/ml-%SCRIPTBASE%.err
# Run in the 'shared' queue (job may share node with other jobs)
#SBATCH -p shared


# Commands to be run:

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
module load gcc/11.2 boost openmpi/4.1.1
mpirun /nobackup/$USER/revbayes/projects/cmake/build-mpi/rb-mpi marginal.Rev

echo "RevBayes job has terminated. Syncing repo:"

# Sync repo with any changes that have happened remotely during the run
git fetch --depth 1 origin

# Commit output files to git
git add %SCRIPTBASE%*.pp
git commit -m "Marginal likelihood output files"
git rebase
git push origin main

# Record usage of tmpDir
du -hs $TMPDIR > /nobackup/$USER/%SCRIPTBASE%/ml-tmpdir_usage.log
