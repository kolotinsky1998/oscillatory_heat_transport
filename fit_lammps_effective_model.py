#!/usr/bin/env python3
"""Fit an effective reservoir ODE to LAMMPS three_bath_observables.dat.

The direct analytic toy model is

    dT/dt = a S T - b (T - Tamb)

where S is the heat-transfer matrix of the nonreciprocal oscillator network.
For a true direct Langevin-reservoir model, a and b can be predicted from the
oscillator friction and the chosen reservoir leakage.  In the atomistic LAMMPS
model, however, a must absorb the finite Ar heat capacity and the Cu-Ar
interface conductance.  This script fits a and b from LAMMPS data.

Usage on the directory containing three_bath_observables.dat:

    python3 fit_lammps_effective_model.py three_bath_observables.dat

Outputs under fitted_effective_model/.
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

from oscillatory_reservoir_python import solve_linear


S = [
    [0.00618868, -0.00777127, -0.00133801],
    [-0.00133801, 0.00618868, -0.00777127],
    [-0.00777127, -0.00133801, 0.00618868],
]

CVV = [
    [0.993811, 0.00777127, 0.00133801],
    [0.00133801, 0.993811, 0.00777127],
    [0.00777127, 0.00133801, 0.993811],
]


def read_observables(path: Path) -> list[list[float]]:
    rows = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            # columns: step time Tamb leakDamp Tar1 Tar2 Tar3 ...
            if len(parts) < 7:
                continue
            rows.append([float(x) for x in parts[:10]])
    return rows


def moving_average(values: list[list[float]], window: int) -> list[list[float]]:
    if window <= 1:
        return [row[:] for row in values]
    half = window // 2
    out = []
    for i in range(len(values)):
        lo = max(0, i - half)
        hi = min(len(values), i + half + 1)
        n = hi - lo
        out.append([sum(values[k][j] for k in range(lo, hi)) / n for j in range(len(values[0]))])
    return out


def matvec(a: list[list[float]], x: list[float]) -> list[float]:
    return [sum(a[i][j] * x[j] for j in range(len(x))) for i in range(len(a))]


def fit_constrained(rows: list[list[float]], window: int = 21) -> tuple[float, float, list[list[float]]]:
    # Extract and smooth [time, Tamb, Tar1, Tar2, Tar3].
    base = [[r[1], r[2], r[4], r[5], r[6]] for r in rows]
    smooth = moving_average(base, window)

    xtx = [[0.0, 0.0], [0.0, 0.0]]
    xty = [0.0, 0.0]
    fit_rows = []

    for k in range(1, len(smooth) - 1):
        t_prev, _, *tp = smooth[k - 1]
        t, tamb, *temp = smooth[k]
        t_next, _, *tn = smooth[k + 1]
        dt = t_next - t_prev
        if dt <= 0:
            continue
        dtemp = [(tn[i] - tp[i]) / dt for i in range(3)]
        st = matvec(S, temp)
        for i in range(3):
            x = [st[i], -(temp[i] - tamb)]
            y = dtemp[i]
            for a in range(2):
                xty[a] += x[a] * y
                for b in range(2):
                    xtx[a][b] += x[a] * x[b]
        fit_rows.append([t, tamb] + temp + dtemp + st)

    coeff = solve_linear(xtx, xty)
    return coeff[0].real, coeff[1].real, fit_rows


def simulate_constrained(
    times: list[float], initial: list[float], tamb: float, a_heat: float, b_leak: float
) -> list[list[float]]:
    def rhs(temp: list[float]) -> list[float]:
        st = matvec(S, temp)
        return [a_heat * st[i] - b_leak * (temp[i] - tamb) for i in range(3)]

    temp = initial[:]
    out = []
    for idx, t in enumerate(times):
        tnp = matvec(CVV, temp)
        out.append([t] + temp[:] + tnp[:])
        if idx == len(times) - 1:
            break
        dt = times[idx + 1] - t
        k1 = rhs(temp)
        k2 = rhs([temp[i] + 0.5 * dt * k1[i] for i in range(3)])
        k3 = rhs([temp[i] + 0.5 * dt * k2[i] for i in range(3)])
        k4 = rhs([temp[i] + dt * k3[i] for i in range(3)])
        temp = [temp[i] + dt * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6.0 for i in range(3)]
    return out


def write_csv(path: Path, header: list[str], rows: list[list[float]]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def write_gnuplot(outdir: Path) -> None:
    (outdir / "plot_fit.gp").write_text(
        """set terminal pngcairo size 1300,850 enhanced font "Arial,12"
set datafile separator comma
set key outside
set grid

set output "fit_reservoir_temperatures.png"
set xlabel "time (ps)"
set ylabel "temperature (K)"
set title "LAMMPS Ar temperatures versus fitted constrained effective model"
plot "lammps_temperatures.csv" using 1:2 with lines lw 1 title "LAMMPS Tar1", \
     "lammps_temperatures.csv" using 1:3 with lines lw 1 title "LAMMPS Tar2", \
     "lammps_temperatures.csv" using 1:4 with lines lw 1 title "LAMMPS Tar3", \
     "fitted_constrained_model.csv" using 1:2 with lines lw 2 title "fit Tar1", \
     "fitted_constrained_model.csv" using 1:3 with lines lw 2 title "fit Tar2", \
     "fitted_constrained_model.csv" using 1:4 with lines lw 2 title "fit Tar3"
"""
    )


def main() -> None:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("three_bath_observables.dat")
    rows = read_observables(path)
    if len(rows) < 10:
        raise SystemExit(f"Not enough rows in {path}")

    outdir = Path("fitted_effective_model")
    outdir.mkdir(exist_ok=True)

    a_heat, b_leak, fit_rows = fit_constrained(rows)
    times = [r[1] for r in rows]
    initial = [rows[0][4], rows[0][5], rows[0][6]]
    tamb = rows[0][2]
    sim = simulate_constrained(times, initial, tamb, a_heat, b_leak)

    lammps_t = [[r[1], r[4], r[5], r[6], r[7], r[8], r[9]] for r in rows]
    write_csv(outdir / "lammps_temperatures.csv", ["time_ps", "Tar1", "Tar2", "Tar3", "Tnp1", "Tnp2", "Tnp3"], lammps_t)
    write_csv(
        outdir / "fitted_constrained_model.csv",
        ["time_ps", "Tar1_fit", "Tar2_fit", "Tar3_fit", "Tnp1_fit", "Tnp2_fit", "Tnp3_fit"],
        sim,
    )
    write_csv(
        outdir / "derivative_fit_samples.csv",
        ["time_ps", "Tamb", "Tar1", "Tar2", "Tar3", "dTar1_dt", "dTar2_dt", "dTar3_dt", "ST1", "ST2", "ST3"],
        fit_rows,
    )
    write_gnuplot(outdir)

    leak_damp_fit = 2.0 / b_leak if b_leak > 0 else float("inf")
    summary = [
        ["a_heat_ps^-1", a_heat],
        ["b_leak_ps^-1", b_leak],
        ["equivalent_temp_leak_time_1_over_b_ps", 1.0 / b_leak if b_leak > 0 else float("inf")],
        ["equivalent_lammps_leakDamp_ps_if_temp_relax_2_over_damp", leak_damp_fit],
    ]
    write_csv(outdir / "fit_summary.csv", ["quantity", "value"], summary)

    print("Constrained effective model fit:")
    print(f"  dT/dt = a_heat*S*T - b_leak*(T-Tamb)")
    print(f"  a_heat = {a_heat:.6g} ps^-1")
    print(f"  b_leak = {b_leak:.6g} ps^-1")
    print(f"  1/b_leak = {1.0 / b_leak if b_leak > 0 else float('inf'):.3g} ps")
    print(f"  equivalent LAMMPS leakDamp ~ 2/b = {leak_damp_fit:.3g} ps")
    print(f"Wrote {outdir}")


if __name__ == "__main__":
    main()

