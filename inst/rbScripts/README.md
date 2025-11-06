## Git checkout on Hamilton
1. Log in to Hamilton using PuTTY
1. Create a personal access token (GitHub user photo → settings → Developer → PAT) that has
   repo read rights for `ms609/neotrans` AND read/write rights for all repos in `neo-trans`
1. Execute `git clone {personal access token}@github.com/ms609/neotrans.git /nobackup/$USER/neotrans/`

## Running RevBayes on Hamilton

1. Log in to Hamilton using PuTTY
1. Navigate to your storage drive: `cd /nobackup/$USER/`
1. Download RevBayes Singularity image, built with MPI
   - Execute `wget https://github.com/revbayes/revbayes/releases/download/v1.2.4/revbayes-v1.2.4-linux64-singularity.simg`
1. Load compiler and OpenMPI modules
   - `module load gcc/11.2 openmpi/4.1.1`
1. Navigate to `cd neotrans`
1. Run the Singularity command through `mpirun`, which will use the local MPI
   - `mpirun -np 2 singularity run -B $TMPDIR:$TMPDIR --app rbmpi ../revbayes-v1.2.4-linux64-singularity.simg rbScripts/script.Rev`
   - `-np 2` is not required when the command is in a SLURM batch script;
     it will be set automatically to match `#SBATCH -n`

## Queueing batch jobs using slurm
1. Create slurm job requests locally in R, using `R/2_RevBayes.R`.
   Check that `options("nt-slurm" = TRUE)` in `R/_setup.R`.
1. Push the files to GitHub, then `cd neotrans`, and `git pull`.
1. Run `sbatch slurm/jobid.sh` to enqueue the job
1. `squeue -u YOUR-USER-ID` will show the status of your currently queued jobs.
