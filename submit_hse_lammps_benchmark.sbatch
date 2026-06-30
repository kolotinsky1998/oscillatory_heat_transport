#!/bin/bash
#SBATCH --job-name=oht-bench
#SBATCH --partition=rocky
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --gpus=1
#SBATCH --constraint="type_a|type_b|type_c"
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err
#SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=your_email@hse.ru

set -euo pipefail

module purge
module load el9
module load lammps/22Jul2025_update3-foss-2025b-kokkos-CUDA-12.9.1

export OMP_NUM_THREADS=1

# Short calibration run.  Make a temporary short input without touching the
# production file.
sed -e 's/variable        equilSteps equal .*/variable        equilSteps equal 10000  # benchmark only/' \
    -e 's/variable        prodSteps equal .*/variable        prodSteps equal 100000  # benchmark only/' \
    in.three_bath_oscillatory_heat.lmp > in.three_bath_oscillatory_heat.bench.lmp

srun --mpi=pmix_v3 lmp -k on g 1 -sf kk \
  -in in.three_bath_oscillatory_heat.bench.lmp
