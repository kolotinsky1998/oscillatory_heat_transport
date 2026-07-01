# Coarse-Grained Bead Open-Chain MD

This setup is an accelerated finite-reservoir demonstration of the PRE2025
heat-redistribution mechanism.  It is not intended as a quantitative
Cu-nanoparticle/Ar model.

## Files

- `make_three_bath_bead_data.py` generates `three_bath_bead.data`.
- `in.three_bath_bead_open_chain.lmp` runs the MD model.
- `plot_bead_open_chain.gp` plots the observables.
- `submit_hse_bead_open_chain.sbatch` is the HSE cluster submission script.

## Model

Each reservoir contains `200` Ar atoms and one heavy bead with mass
`8000 amu`.  The bead mass ratio is

```text
M_bead / m_Ar = 8000 / 39.948 ~= 200
```

which keeps the bead in a reasonable heavy-particle Langevin regime while
being much faster than a full Cu nanoparticle.

The production run has no thermostat.  Langevin thermostats are used only in
the initial equilibration stage and are removed before `reset_timestep 0`.

The open-chain force matrix is

```text
A = [[1.0, 0.6, 0.0],
     [0.3, 1.0, 0.6],
     [0.0, 0.3, 1.0]]
```

with `kforce = kperp = 0.20 eV/A^2`.

This matrix is conservative in the PRE2025 sense:

```text
d_i A_ij = d_j A_ji,  d = [1, 2, 4]
```

so the expected stationary finite-bath ratio is

```text
T1:T2:T3 = 1 : 1/2 : 1/4
```

or, normalized to mean `140 K`,

```text
T1 ~= 240 K, T2 ~= 120 K, T3 ~= 60 K
```

The actual MD values can differ because the controller may do net work and the
Ar reservoirs are deliberately tiny.  Use `bead_open_chain_total_energy.png`
and `bead_open_chain_feedback_power.png` to diagnose that.

## Run

```bash
python3 make_three_bath_bead_data.py --output three_bath_bead.data
lmp -in in.three_bath_bead_open_chain.lmp
gnuplot plot_bead_open_chain.gp
```

On HSE:

```bash
sbatch submit_hse_bead_open_chain.sbatch
```

Main output:

```text
bead_open_chain_observables.dat
```

Key plots:

- `bead_open_chain_reservoir_temperatures.png`
- `bead_open_chain_bead_temperatures.png`
- `bead_open_chain_total_energy.png`
- `bead_open_chain_feedback_power.png`
