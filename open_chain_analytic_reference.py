#!/usr/bin/env python3
"""Analytic reference for the open-chain conservative nonreciprocal test."""

from __future__ import annotations

from pathlib import Path

import oscillatory_reservoir_python as model


F_OPEN = [
    [1.00, 0.14, 0.00],
    [0.06, 1.00, 0.14],
    [0.00, 0.06, 1.00],
]


def main() -> None:
    old_f = [row[:] for row in model.F]
    old_gamma = model.GAMMA
    model.F[:] = [row[:] for row in F_OPEN]
    model.GAMMA = 1.10578
    try:
        s, c_vv = model.heat_transfer_matrix()
        eigs = model.eig3(s)
    finally:
        model.F[:] = old_f
        model.GAMMA = old_gamma

    # Conservativity condition d_i F_ij = d_j F_ji.
    # Choose d1=1.  Then d2=d1*f12/f21 and d3=d2*f23/f32.
    d1 = 1.0
    d2 = d1 * F_OPEN[0][1] / F_OPEN[1][0]
    d3 = d2 * F_OPEN[1][2] / F_OPEN[2][1]
    d = [d1, d2, d3]
    db = [1.0 / x for x in d]
    mean_temp = 140.0
    scale = 3.0 * mean_temp / sum(db)
    db_scaled = [scale * x for x in db]

    text = []
    text.append("# Open-Chain Analytic Reference\n")
    text.append("Open-chain force matrix:\n")
    text.append("```text")
    for row in F_OPEN:
        text.append(" ".join(f"{x: .6g}" for x in row))
    text.append("```\n")
    text.append("Symmetrizing diagonal `D = diag(d_i)` from `d_i F_ij = d_j F_ji`:\n")
    text.append(f"```text\nd = {d}\n```\n")
    text.append("Detailed balance in the original variables requires `D.T = const`, so `T_i proportional 1/d_i`:\n")
    text.append(f"```text\nT ratio = {db}\nnormalized to mean 140 K: {db_scaled}\n```\n")
    text.append("Heat-transfer matrix S for gamma*=1.10578:\n")
    text.append("```text")
    for row in s:
        text.append(" ".join(f"{x: .8g}" for x in row))
    text.append("```\n")
    text.append("Eigenvalues of S:\n")
    text.append("```text")
    for z in eigs:
        text.append(f"{z.real:.8g} {z.imag:+.8g} i")
    text.append("```\n")
    text.append("Cvv = I - S, so analytic Cu kinetic temperatures are Tnp = Cvv . Tar:\n")
    text.append("```text")
    for row in c_vv:
        text.append(" ".join(f"{x: .8g}" for x in row))
    text.append("```\n")

    out = Path("open_chain_analytic_reference.md")
    out.write_text("\n".join(text))
    print(out.read_text())


if __name__ == "__main__":
    main()
