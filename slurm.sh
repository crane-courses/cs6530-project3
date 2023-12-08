#!/bin/bash
#
#SBATCH --account=soc-np
#SBATCH --partition=soc-np
#
#SBATCH --job-name=b-skip-list
#SBATCH -o /uufs/chpc.utah.edu/common/home/sullivan-group1/logs/slurm/b-skip-list/b-skip-list-%j-%N.out
#SBATCH -e /uufs/chpc.utah.edu/common/home/sullivan-group1/logs/slurm/b-skip-list/b-skip-list-%j-%N.err
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#
#SBATCH --mail-user=alex.crane@utah.edu
#SBATCH --mail-type=ALL
#
#SBATCH --time=7-00:00:00

set -e

# set paths
SHARED_DIR="/uufs/chpc.utah.edu/common/home/sullivan-group1"
OUTPUT_DIR="${SHARED_DIR}/b-skip-list/results"

# assumes that you already have the repo cloned and compiled
WORK_DIR="${HOME}/dev/b-skip-list"
cd "${WORK_DIR}"

module load gcc/11.2.0
# make clean
# make
./basic 10000

# collect Slurm stats
echo "===="
echo "Slurm environment variables:"
env | grep SLURM
echo "----"
scontrol show job ${SLURM_JOB_ID}     # Print out final statistics about resource uses before job exits
echo "----"
sstat -p --format 'JobID,MaxRSS,AveCPU,MaxDiskRead,MaxDiskWrite' -j ${SLURM_JOB_ID}.batch
echo "===="

exit 0