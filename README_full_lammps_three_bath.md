# Full LAMMPS three-bath oscillatory heat-transport model

This is the fuller atomistic/coarse-atom MD version inspired by Loos et al.,
Scientific Reports 13, 4517 (2023).  It uses explicit Ar reservoirs, Cu
nanoparticles, fixed Cu walls, LJ interactions, and an external nonreciprocal
COM feedback controller.

## Files

- `make_three_bath_data.py` - generates the LAMMPS data file.
- `three_bath_cu_ar.data` - generated initial configuration.
- `in.three_bath_oscillatory_heat.lmp` - main LAMMPS input.
- `plot_three_bath_observables.gp` - gnuplot script for result figures.

The older `in.oscillatory_reservoirs.lmp` is only a three-particle toy model.
Use `in.three_bath_oscillatory_heat.lmp` for the many-core MD run.

## Generate Data

```bash
python3 make_three_bath_data.py --output three_bath_cu_ar.data
```

With default settings this creates about `65k` atoms:

- three Ar reservoirs, about `12812` Ar atoms each;
- three Cu nanoparticles, about `1055` Cu atoms each;
- fixed Cu walls at the box boundaries and between reservoirs.

The default nanoparticle radius is `14.4 A`, matching the radius scale quoted
in the paper.  For a faster test, reduce `--compartment-length`, `--box-y`,
`--box-z`, or `--np-radius`.

## Run

Example MPI run:

```bash
mpirun -np 16 lmp -in in.three_bath_oscillatory_heat.lmp
```

The input uses:

- `units metal`;
- timestep `0.001 ps = 1 fs`;
- LJ parameters from the paper:
  - Ar-Ar: `epsilon = 0.0104 eV`, `sigma = 3.405 A`;
  - Cu-Cu: `epsilon = 0.4093 eV`, `sigma = 2.338 A`;
  - Ar-Cu: Lorentz-Berthelot mixed values.
- equilibration Langevin damping `0.1 ps`;
- production leakage damping `200 ps` by default;
- equilibration `50 ps`;
- production `15 ns` by default.

The current production setting is already:

```lammps
variable        prodSteps equal 15000000
```

which gives `15 ns` at `1 fs` timestep.

## Nonreciprocal Controller

The controller acts on the nanoparticle centers of mass.  The total COM force is
distributed uniformly over atoms in each nanoparticle:

```text
Fcom_i = -kforce sum_j F_ij (X_j - X_j0)
```

with

```text
F = [[1.00, 0.14, 0.06],
     [0.06, 1.00, 0.14],
     [0.14, 0.06, 1.00]]
```

and `kforce = 0.01 eV/A^2`, corresponding to `1 eV/nm^2`.

## Reservoir Temperatures

During production, the Ar reservoir temperatures are not prescribed by an
oscillatory target.  The measured reservoir temperatures are `Tar1`, `Tar2`,
and `Tar3`, and they evolve due to heat exchange with the Cu nanoparticles.

The only thermostat in production is a weak external leakage bath:

```text
fix leak1 ar1 langevin Tamb Tamb leakDamp ...
fix leak2 ar2 langevin Tamb Tamb leakDamp ...
fix leak3 ar3 langevin Tamb Tamb leakDamp ...
```

This represents the additional heat outflow term that prevents unlimited
temperature growth:

```text
dT/dt = heat_exchange_with_nanoparticles - leakage_rate*(T - Tamb)
```

The default leakage parameters are:

```text
Tamb = 110 K
leakDamp = 200 ps
```

Increase `leakDamp` for weaker leakage; decrease it for stronger stabilization.
If `leakDamp` is too small, it will simply pin the Ar reservoirs to `Tamb` and
wash out the oscillations.

The initial reservoir temperatures are imposed only before production:

```text
Tinit1 = 173.5 K
Tinit2 = 143.5 K
Tinit3 = 125.5 K
```

## Main Output

LAMMPS writes:

```text
three_bath_observables.dat
```

Columns:

```text
step time_ps Tamb leakDamp
Tar1 Tar2 Tar3 Tnp1 Tnp2 Tnp3
X1 X2 X3 u1 u2 u3
Fcom1 Fcom2 Fcom3
Qleak1 Qleak2 Qleak3
```

Plot:

```bash
gnuplot plot_three_bath_observables.gp
```

Expected figures:

- `three_bath_temperatures.png` - main self-consistent reservoir-temperature result;
- `three_bath_nanoparticle_temperatures.png` - Cu nanoparticle kinetic temperatures;
- `three_bath_com_displacements.png` - nanoparticle COM displacement dynamics;
- `three_bath_feedback_forces.png` - nonreciprocal controller forces;
- `three_bath_leakage_energy.png` - cumulative energy exchanged with the weak external leakage baths.
