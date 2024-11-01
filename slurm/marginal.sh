#!/bin/bash
 
# Request resources:
#SBATCH -n 16          # number of tasks, up to 128. ## Replacing "SBATCH -c 16" : # number of CPU cores, one per thread, up to 128
#SBATCH --mem=64G      # memory required, up to 250G on standard nodes
#SBATCH --time=2-23:45:0  # time limit for job (format:  days-hours:minutes:seconds)
#SBATCH --gres=tmp:128G # temporary disk space required on the compute node ($TMPDIR), up to 400G
#SBATCH --job-name=ml-%SCRIPTBASE%
#SBATCH --output=ml-%SCRIPTBASE%.out
#SBATCH --error=ml-%SCRIPTBASE%.err
# Run in the 'shared' queue (job may share node with other jobs)
#SBATCH -p shared

# Commands to be run:

# Set git identifier
git config --global user.name "RevBayes analysis"

# Check out analysis repo
git clone --depth 1 https://${MKNT_WRITE}@github.com/mk-nt/%SCRIPTBASE%.git /nobackup/$USER/%SCRIPTBASE%
cd /nobackup/$USER/%SCRIPTBASE%
git checkout main
git pull origin main --rebase

# Run RevBayes
module load gcc/11.2 boost openmpi/4.1.1
mpirun ../revbayes/projects/cmake/build-mpi/rb-mpi marginal.Rev

echo "RevBayes job has terminated. Syncing repo:"

# Sync repo with any changes that have happened remotely during the run
git fetch --depth 1 origin
git rebase origin

# Commit output files to git
git add %SCRIPTBASE%*.pp
git commit -m "Marginal likelihood output files"
git push origin main
