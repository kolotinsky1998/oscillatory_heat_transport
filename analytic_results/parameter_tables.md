# Analytic LAMMPS Mapping Tables

## Physical Parameters

| quantity | value | unit |
|---|---:|---|
| N_Cu per nanoparticle | 1055 | atoms |
| M_NP | 1.11324e-22 | kg |
| R_NP | 1.44 | nm |
| Ar number density | 29.765 | nm^-3 |
| Ar mass density | 1974.47 | kg/m^3 |
| K diagonal | 0.01 | eV/A^2 |
| K diagonal | 0.160218 | N/m |
| t0=sqrt(M/K) | 26.3597 | ps |
| Gamma_ref | 4.67e-12 | kg/s |
| Gamma_star | 1.10578 | dimensionless |
| M/Gamma_ref | 23.8382 | ps |
| M/(2 Gamma_ref) | 11.9191 | ps |
| leakDamp | 200 | ps |
| r_star=2 t0/leakDamp | 0.263597 | dimensionless |

## Epstein Drag Estimates

| T (K) | Gamma delta=1 (kg/s) | M/Gamma (ps) | Gamma delta=1.39 (kg/s) | M/Gamma (ps) |
|---:|---:|---:|---:|---:|
| 110 | 4.1409e-12 | 26.88 | 5.7559e-12 | 19.34 |
| 140 | 4.6716e-12 | 23.83 | 6.4935e-12 | 17.14 |
| 170 | 5.1479e-12 | 21.63 | 7.1555e-12 | 15.56 |

## Heat Transfer Matrix S for Current Physical Mapping

Rows define `J = gamma S T` in the dimensionless oscillator model.

| row | col1 | col2 | col3 |
|---:|---:|---:|---:|
| 1 | 0.00618868 | -0.00777127 | -0.00133801 |
| 2 | -0.00133801 | 0.00618868 | -0.00777127 |
| 3 | -0.00777127 | -0.00133801 | 0.00618868 |

## Effective Cu Kinetic Temperature Response Cvv = I - S

Analytic nanoparticle kinetic temperatures are `Tnp = Cvv . Tar`.

| row | col1 | col2 | col3 |
|---:|---:|---:|---:|
| 1 | 0.993811 | 0.00777127 | 0.00133801 |
| 2 | 0.00133801 | 0.993811 | 0.00777127 |
| 3 | 0.00777127 | 0.00133801 | 0.993811 |

## Eigenvalues

| matrix | eigenvalue |
|---|---:|
| S | -0.0029206 -1.25354e-32 i |
| S | 0.0107433 +0.00557137 i |
| S | 0.0107433 -0.00557137 i |
| reservoir A | -0.328187 +1.31372e-46 i, decay=80.3 ps |
| reservoir A | -0.0260027 +0.123214 i, period=1344.2 ps, decay=1013.7 ps |
| reservoir A | -0.0260027 -0.123214 i, period=1344.2 ps, decay=1013.7 ps |

## Heat-Rate Scale alphaHeat

The toy analytic value `kappa=20` is kept, but the atomistic reservoir has a large finite heat capacity.
Use `alphaHeat` as a calibration multiplier in

`dT/dtau = alphaHeat*kappa*GammaStar*S*T - rStar*(T-Tamb)`.

| alphaHeat | final mean T at 1800 ps for current LAMMPS initial condition (K) |
|---:|---:|
| 1 | 88.351 |
| 0.3 | 102.468 |
| 0.1 | 107.369 |
| 0.03 | 109.197 |
| 0.01 | 109.731 |
