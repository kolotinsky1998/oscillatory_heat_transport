#!/usr/bin/env python3
"""Estimate heat-flux intensity and reservoir relaxation time for the open chain.

The calculation follows PRE2025 Eq. (30) and Eq. (36):

    S_ij = delta_ij - gamma * int omega^2 |G_ij(omega)|^2 d omega / pi
    dT/dt = kappa_heat * gamma * S T

For the atomistic Ar reservoirs, kappa_heat is obtained from the finite heat
capacity of one reservoir.  The output is intended as an order-of-magnitude
guide for how long the LAMMPS open-chain run must be.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np


KB = 1.380649e-23
AMU = 1.66053906660e-27
EV_TO_J = 1.602176634e-19
EV_PER_A2_TO_N_PER_M = 16.0217662

F_OPEN = np.array(
    [
        [1.00, 0.14, 0.00],
        [0.06, 1.00, 0.14],
        [0.00, 0.06, 1.00],
    ],
    dtype=float,
)

N_AR = 12812
N_CU_NP = 1055
M_CU = 63.546 * AMU
M_NP = N_CU_NP * M_CU

KFORCE_EV_A2 = 0.01
K_SI = KFORCE_EV_A2 * EV_PER_A2_TO_N_PER_M

# Epstein estimate used in analytic_results/parameter_tables.md near 140 K.
GAMMA_SI = 4.67e-12

# Characteristic oscillator units used in PRE2025-style dimensionless model.
T0_S = math.sqrt(M_NP / K_SI)
GAMMA_STAR = GAMMA_SI * T0_S / M_NP

# LAMMPS reservoirs are 3D monatomic Ar fluids.  For a dense LJ fluid the
# constant-volume heat capacity is not exactly 3/2 NkB, so we report both the
# kinetic lower estimate and a conservative 3 NkB estimate.
CV_KINETIC = 1.5 * N_AR * KB
CV_FLUID_EST = 3.0 * N_AR * KB


def heat_transfer_matrix_lyapunov(gamma_star: float) -> np.ndarray:
    """Compute S from the equivalent stationary covariance problem.

    Directly integrating Eq. (30) over a finite frequency interval converges
    slowly because the diagonal terms have a 1/omega^2 tail.  Solving the
    Lyapunov equation gives the same Cvv response without tail truncation.
    """
    n = F_OPEN.shape[0]
    dim = 2 * n
    a = np.zeros((dim, dim), dtype=float)
    a[:n, n:] = np.eye(n)
    a[n:, :n] = -F_OPEN
    a[n:, n:] = -gamma_star * np.eye(n)

    c_vv = np.zeros((n, n), dtype=float)
    lhs = np.kron(np.eye(dim), a) + np.kron(a, np.eye(dim))
    for j in range(n):
        q = np.zeros((dim, dim), dtype=float)
        q[n + j, n + j] = 2.0 * gamma_star
        cov = np.linalg.solve(lhs, -q.reshape(dim * dim)).reshape((dim, dim))
        c_vv[:, j] = np.diag(cov[n:, n:])
    return np.eye(n) - c_vv


def stationary_temperature_ratio() -> np.ndarray:
    # Conservativity condition d_i F_ij = d_j F_ji.
    d = np.array([1.0, F_OPEN[0, 1] / F_OPEN[1, 0], (F_OPEN[0, 1] / F_OPEN[1, 0]) * F_OPEN[1, 2] / F_OPEN[2, 1]])
    ratio = 1.0 / d
    return ratio / ratio.mean()


def format_matrix(a: np.ndarray) -> str:
    return "\n".join("  " + " ".join(f"{x: .8e}" for x in row) for row in a)


def main() -> None:
    s = heat_transfer_matrix_lyapunov(GAMMA_STAR)
    eig = np.linalg.eigvals(s)
    ratio = stationary_temperature_ratio()

    # Physical heat power scale for Langevin friction:
    #   J_i = (Gamma/M) kB (S T)_i.
    # The reservoir loses this heat, so dT_res/dt = -(Gamma/M) kB S T / C_V.
    pref_kin = (GAMMA_SI / M_NP) * KB / CV_KINETIC
    pref_fluid = (GAMMA_SI / M_NP) * KB / CV_FLUID_EST

    rates_kin = pref_kin * eig
    rates_fluid = pref_fluid * eig

    nonzero = [z for z in eig if abs(z) > 1e-8]
    # Representative heat flux per reservoir for a 1 K deviation along the
    # slowest nonzero mode.  Since J = gamma kB S T, eigenvalue magnitude gives
    # the scale.
    eig_nonzero_sorted = sorted(nonzero, key=lambda z: abs(z.real))
    lambda_slow = eig_nonzero_sorted[0].real
    heat_flux_per_k = abs((GAMMA_SI / M_NP) * KB * lambda_slow)
    temp_rate_per_k_fluid = heat_flux_per_k / CV_FLUID_EST
    slow_s = 1.0 / (pref_fluid * abs(lambda_slow))

    out = []
    out.append("# Open-chain heat-flux estimate from PRE2025\n")
    out.append("Force matrix A:")
    out.append("```text")
    out.append(format_matrix(F_OPEN))
    out.append("```\n")
    out.append(f"N_Ar per reservoir: {N_AR}")
    out.append(f"N_Cu per nanoparticle: {N_CU_NP}")
    out.append(f"M_NP: {M_NP:.8e} kg")
    out.append(f"k: {KFORCE_EV_A2:g} eV/A^2 = {K_SI:.8e} N/m")
    out.append(f"t0 = sqrt(M/k): {T0_S * 1e12:.6g} ps")
    out.append(f"gamma: {GAMMA_SI:.8e} kg/s")
    out.append(f"gamma*: {GAMMA_STAR:.8g}\n")
    out.append("Heat-transfer matrix S from PRE2025 Eq. (30), computed via equivalent covariance problem:")
    out.append("```text")
    out.append(format_matrix(s))
    out.append("```\n")
    out.append("Eigenvalues of S:")
    out.append("```text")
    for z in eig:
        out.append(f"{z.real:.10e} {z.imag:+.10e} i")
    out.append("```\n")
    out.append("Stationary open-chain temperature ratio, normalized to mean 1:")
    out.append("```text")
    out.append(" ".join(f"{x:.8g}" for x in ratio))
    out.append("```\n")
    out.append("Finite-reservoir conversion:")
    out.append(f"- Cv kinetic only = 3/2 N kB = {CV_KINETIC:.8e} J/K")
    out.append(f"- Cv fluid estimate = 3 N kB = {CV_FLUID_EST:.8e} J/K")
    out.append(f"- (gamma/M) kB / Cv kinetic = {pref_kin:.8e} 1/s")
    out.append(f"- (gamma/M) kB / Cv fluid = {pref_fluid:.8e} 1/s\n")
    out.append("Reservoir relaxation rates A = -(gamma/M)(kB/Cv) S:")
    out.append("```text")
    out.append("kinetic Cv:")
    for z in rates_kin:
        rate = -z
        tau = math.inf if abs(rate.real) < 1e-30 else -1.0 / rate.real
        out.append(f"  rate={rate.real:.8e} {rate.imag:+.8e} i 1/s, tau={tau * 1e9:.6g} ns")
    out.append("fluid Cv:")
    for z in rates_fluid:
        rate = -z
        tau = math.inf if abs(rate.real) < 1e-30 else -1.0 / rate.real
        out.append(f"  rate={rate.real:.8e} {rate.imag:+.8e} i 1/s, tau={tau * 1e9:.6g} ns")
    out.append("```\n")
    out.append(f"Slow heat-flux scale per 1 K along the slow mode: {heat_flux_per_k:.8e} W")
    out.append(f"Corresponding fluid-Cv temperature rate per 1 K: {temp_rate_per_k_fluid:.8e} K/s")
    out.append(f"Characteristic slow relaxation time with fluid Cv: {slow_s * 1e9:.6g} ns")
    out.append(f"Rule of thumb: 3 tau = {3 * slow_s * 1e9:.6g} ns, 5 tau = {5 * slow_s * 1e9:.6g} ns")

    path = Path("open_chain_heat_flux_estimate.md")
    path.write_text("\n".join(out))
    print(path.read_text())


if __name__ == "__main__":
    main()
