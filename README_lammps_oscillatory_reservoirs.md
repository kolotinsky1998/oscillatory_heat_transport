# LAMMPS oscillatory-reservoir demo

## What is taken from the MD paper

The paper `s41598-023-31583-y (1).pdf` uses LAMMPS all-atom MD for two copper nanoparticles in separate argon baths. The key technical ingredients are:

- Cu nanoparticle radius about `1.4 nm`, made of `186` Cu atoms.
- Two independent Ar reservoirs with `25280` Ar atoms each.
- Reservoir temperatures `TC = 100 K`, `TH = 120 K`.
- Static harmonic traps with stiffness `kappa = 1 eV/nm^2`.
- External feedback/nonreciprocal COM forces, e.g. `kappaC (XH - XC)` and `kappaH (XC - XH)`.
- NVT reservoirs with Nose-Hoover coupling time `1 ps`; Langevin thermostats were also checked in the paper.
- MD timestep `1 fs`, sampled every `1 ps`, total trajectory `15-30 ns`.
- Effective COM model:
  `m X'' + gamma X' = trap force + nonreciprocal force + noise`,
  with fitted `gamma approximately 3.5e-12 kg/s`.

The LAMMPS file here keeps the effective underdamped Langevin/feedback layer, not the full Cu/Ar atomistic bath. That is the practical version for a one-core demonstration of the temperature oscillation mechanism.

## Files

- `in.oscillatory_reservoirs.lmp` - LAMMPS input.
- `plot_reservoir_temperatures.gp` - gnuplot script for PNG plots.

## Run

```bash
lmp -in in.oscillatory_reservoirs.lmp
gnuplot plot_reservoir_temperatures.gp
```

or explicitly on one MPI rank:

```bash
mpirun -np 1 lmp -in in.oscillatory_reservoirs.lmp
```

## Output

The main data file is:

```text
reservoir_temperatures.dat
```

Columns are:

```text
step time T1 T2 T3 x1 x2 x3 F1 F2 F3 Qlangevin1 Qlangevin2 Qlangevin3
```

The first plot to inspect is:

```text
reservoir_temperatures.png
```

It should show damped oscillatory reservoir temperatures. The damping comes from the added heat outflow term in the reservoir equation.

