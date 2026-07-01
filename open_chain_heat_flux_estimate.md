# Open-chain heat-flux estimate from PRE2025

Force matrix A:
```text
   1.00000000e+00  1.40000000e-01  0.00000000e+00
   6.00000000e-02  1.00000000e+00  1.40000000e-01
   0.00000000e+00  6.00000000e-02  1.00000000e+00
```

N_Ar per reservoir: 12812
N_Cu per nanoparticle: 1055
M_NP: 1.11324249e-22 kg
k: 0.01 eV/A^2 = 1.60217662e-01 N/m
t0 = sqrt(M/k): 26.3597 ps
gamma: 4.67000000e-12 kg/s
gamma*: 1.1057751

Heat-transfer matrix S from PRE2025 Eq. (30), computed via equivalent covariance problem:
```text
   3.41292548e-03 -7.90616815e-03 -1.33757512e-04
  -1.45215333e-03  6.77671555e-03 -7.90616815e-03
  -4.51243585e-06 -1.45215333e-03  3.41292548e-03
```

Eigenvalues of S:
```text
1.0165073331e-02 +0.0000000000e+00 i
3.4374931897e-03 +0.0000000000e+00 i
-6.5813315097e-17 +0.0000000000e+00 i
```

Stationary open-chain temperature ratio, normalized to mean 1:
```text
1.8607595 0.79746835 0.34177215
```

Finite-reservoir conversion:
- Cv kinetic only = 3/2 N kB = 2.65333125e-19 J/K
- Cv fluid estimate = 3 N kB = 5.30666250e-19 J/K
- (gamma/M) kB / Cv kinetic = 2.18282498e+06 1/s
- (gamma/M) kB / Cv fluid = 1.09141249e+06 1/s

Reservoir relaxation rates A = -(gamma/M)(kB/Cv) S:
```text
kinetic Cv:
  rate=-2.21885760e+04 +0.00000000e+00 i 1/s, tau=45068.2 ns
  rate=-7.50344602e+03 +0.00000000e+00 i 1/s, tau=133272 ns
  rate=1.43658948e-10 +0.00000000e+00 i 1/s, tau=-6.96093e+18 ns
fluid Cv:
  rate=-1.10942880e+04 +0.00000000e+00 i 1/s, tau=90136.5 ns
  rate=-3.75172301e+03 +0.00000000e+00 i 1/s, tau=266544 ns
  rate=7.18294742e-11 +0.00000000e+00 i 1/s, tau=-1.39219e+19 ns
```

Slow heat-flux scale per 1 K along the slow mode: 1.99091278e-15 W
Corresponding fluid-Cv temperature rate per 1 K: 3.75172301e+03 K/s
Characteristic slow relaxation time with fluid Cv: 266544 ns
Rule of thumb: 3 tau = 799633 ns, 5 tau = 1.33272e+06 ns