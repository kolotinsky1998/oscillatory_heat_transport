#!/bin/bash
#SBATCH --job-name=oht-bead
#SBATCH --partition=rocky
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --constraint="type_a|type_b|type_c"
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail

module purge
module load el9
module load lammps/22Jul2025_update3-foss-2025b-kokkos

export OMP_NUM_THREADS=1

python3 make_three_bath_bead_data.py --output three_bath_bead.data

srun -n ${SLURM_NTASKS} lmp -in in.three_bath_bead_open_chain.lmp

gnuplot plot_bead_open_chain.gp
