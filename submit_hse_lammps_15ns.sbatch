#!/bin/bash
#SBATCH --job-name=oht-15ns
#SBATCH --partition=rocky
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --gpus=1
#SBATCH --constraint="type_a|type_b|type_c"
#SBATCH --time=36:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err
#SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=your_email@hse.ru

set -euo pipefail

module purge
module load el9
module load lammps/22Jul2025_update3-foss-2025b-kokkos-CUDA-12.9.1

export OMP_NUM_THREADS=1

# Production: 15 ns at timestep 1 fs.
srun --mpi=pmix_v3 lmp -k on g 1 -sf kk \
  -in in.three_bath_oscillatory_heat.lmp
