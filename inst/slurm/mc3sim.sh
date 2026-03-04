#!/bin/bash
 
# Request resources:
#SBATCH -n 16          # number of tasks, up to 128. ## Replacing "SBATCH -c 16" : # number of CPU cores, one per thread, up to 128
#SBATCH --mem=24M      # memory required, up to 250G on standard nodes
#SBATCH --time=23:45:0  # time limit for job (format:  days-hours:minutes:seconds)
#SBATCH --gres=tmp:128G # temporary disk space required on the compute node ($TMPDIR), up to 400G
#SBATCH --job-name=%SIMNAME%_%SCRIPTID%
#SBATCH --output=/nobackup/%u/neotrans/%SIMNAME%_%SCRIPTID%.out
#SBATCH --error=/nobackup/%u/neotrans/%SIMNAME%_%SCRIPTID%.err
#SBATCH -p shared
#SBATCH --export=ALL

# Set git identifier
git config --global user.name "RevBayes simulation analysis"

# Check out analysis repo - clone if doesn't exist, pull if it does
cd /nobackup/$USER/neotrans/
git pull origin main --rebase

# Run RevBayes
module load gcc/11.2 boost openmpi/4.1.1
start_time=$(date +%s)
mpirun /nobackup/$USER/revbayes/projects/cmake/build-mpi/rb-mpi inst/rbScripts/sim-mc3.Rev inst/%SIMDIR1%/%SIMDIR2% %SCRIPTID% 333
end_time=$(date +%s)
echo "RevBayes job terminated at $end_time."


# Sync repo with any changes that have happened remotely during the run
git fetch --depth 1 origin

cd /nobackup/$USER/neotrans/inst/%SIMDIR1%/%SIMDIR2%
if ls *.trees 1> /dev/null 2>&1; then
  # Record usage of tmpDir
  du -hs $TMPDIR > /nobackup/$USER/neotrans/inst/%SIMDIR1%/%SIMDIR2%/mc3-tmpdir_usage.log

  for file in *.trees; do
    tar -czf "${file%.trees}.tar.gz" "$file"
  done
  
  # Commit output files to git
  git add ./%SCRIPTID%_run_*.tar.gz
  git commit -m "MCMCMC trees: %SIMDIR1% %SIMDIR2%: %SCRIPTID%"
  git rebase main
  git push origin main
else 
  echo "No .trees files found. No output to commit."
fi
