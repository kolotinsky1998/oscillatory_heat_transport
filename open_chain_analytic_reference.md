# Open-Chain Analytic Reference

Open-chain force matrix:

```text
 1  0.14  0
 0.06  1  0.14
 0  0.06  1
```

Symmetrizing diagonal `D = diag(d_i)` from `d_i F_ij = d_j F_ji`:

```text
d = [1.0, 2.3333333333333335, 5.4444444444444455]
```

Detailed balance in the original variables requires `D.T = const`, so `T_i proportional 1/d_i`:

```text
T ratio = [1.0, 0.42857142857142855, 0.18367346938775506]
normalized to mean 140 K: [260.50632911392404, 111.64556962025316, 47.84810126582277]
```

Heat-transfer matrix S for gamma*=1.10578:

```text
 0.0034128954 -0.0079060987 -0.0001337555
-0.0014521406  0.0067766561 -0.0079060987
-4.5123678e-06 -0.0014521406  0.0034128954
```

Eigenvalues of S:

```text
0.010164984 -4.0453136e-36 i
0.0017187313 +0.0022171116 i
0.0017187313 -0.0022171116 i
```

Cvv = I - S, so analytic Cu kinetic temperatures are Tnp = Cvv . Tar:

```text
 0.9965871  0.0079060987  0.0001337555
 0.0014521406  0.99322334  0.0079060987
 4.5123678e-06  0.0014521406  0.9965871
```
