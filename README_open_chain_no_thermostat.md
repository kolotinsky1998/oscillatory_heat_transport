# Open-chain no-thermostat control test

This is a control LAMMPS run for the conservative nonreciprocal force matrix

```text
F = [[1.00, 0.14, 0.00],
     [0.06, 1.00, 0.14],
     [0.00, 0.06, 1.00]]
```

The production run has no Ar thermostats.  The Ar reservoirs are equilibrated
to the same initial temperature, then all Ar thermostats are disabled.  Any
redistribution of reservoir temperatures must come from heat exchange through
the Cu nanoparticles and the open-chain controller.

## Files

- `in.three_bath_open_chain_no_thermostat.lmp` - LAMMPS input.
- `submit_hse_open_chain_no_thermostat.sbatch` - HSE HPC GPU/KOKKOS SLURM script.
- `plot_open_chain_no_thermostat.gp` - plots for the output.
- `open_chain_analytic_reference.py` - analytic detailed-balance reference.
- `open_chain_analytic_reference.md` - generated reference tables.

## Run

```bash
sbatch submit_hse_open_chain_no_thermostat.sbatch
```

Default production length:

```text
prodSteps = 2000000
```

which is `2 ns` at `1 fs`.

## Output

Main output:

```text
open_chain_no_thermostat_observables.dat
```

Columns:

```text
step time_ps
Tar1 Tar2 Tar3
Tnp1 Tnp2 Tnp3
X1 X2 X3
u1 u2 u3
Fcom1 Fcom2 Fcom3
```

Plot:

```bash
gnuplot plot_open_chain_no_thermostat.gp
```

Expected figures:

- `open_chain_reservoir_temperatures.png`
- `open_chain_nanoparticle_temperatures.png`
- `open_chain_com_displacements.png`
- `open_chain_feedback_forces.png`

## Analytic stationary ratio

The matrix is conservative because there exists a positive diagonal matrix
`D = diag(d_i)` such that `D F` is symmetric:

```text
d = [1, 2.3333333333, 5.4444444444]
```

Detailed balance in the original variables requires

```text
D T = const
```

so the expected temperature ratio is

```text
T_i proportional 1/d_i = [1, 0.4285714286, 0.1836734694]
```

Normalized to mean `140 K`, this gives

```text
[260.5 K, 111.6 K, 47.85 K]
```

The exact absolute final temperatures in finite atomistic MD will depend on
the total conserved energy and on how much energy is stored in Cu internal and
COM modes, but this ratio is the reference direction for the stationary
nonuniform distribution.

