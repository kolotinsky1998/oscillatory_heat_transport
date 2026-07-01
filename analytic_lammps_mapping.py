#!/usr/bin/env python3
"""Analytic mapping between the Cu/Ar LAMMPS system and the 3-oscillator model.

Outputs:
  analytic_results/parameter_tables.md
  analytic_results/epstein_gamma.csv
  analytic_results/mapping_summary.csv
  analytic_results/current_lammps_analytic.csv
  analytic_results/current_lammps_reservoirs.svg
  analytic_results/current_lammps_nanoparticles.svg
  analytic_results/centered_excitation_analytic.csv
  analytic_results/centered_excitation_reservoirs.svg
  analytic_results/centered_excitation_nanoparticles.svg
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

from oscillatory_reservoir_python import F, N, eig3, heat_transfer_matrix


OUT = Path("analytic_results")

KB = 1.380649e-23
AMU = 1.66053906660e-27
EV_PER_A2_TO_N_PER_M = 16.0217662

M_AR = 39.948 * AMU
M_CU = 63.546 * AMU
N_CU = 1055
M_NP = N_CU * M_CU
R_NP = 14.4e-10

N_AR = 12812
AR_BOX_SIDE_A = 95.0 - 2.0 * 9.0
AR_EXCLUSION_RADIUS_A = 14.4 + 4.0
V_AR = (AR_BOX_SIDE_A**3 - 4.0 * math.pi * AR_EXCLUSION_RADIUS_A**3 / 3.0) * 1e-30
RHO_AR = N_AR * M_AR / V_AR
N_DENSITY_AR_NM3 = N_AR / (V_AR * 1e27)

KFORCE_EV_A2 = 0.01
K_DIAG = KFORCE_EV_A2 * EV_PER_A2_TO_N_PER_M
T0 = math.sqrt(M_NP / K_DIAG)

GAMMA_REF = 4.67e-12  # Epstein estimate near 140 K, delta about 1
GAMMA_STAR = (GAMMA_REF / M_NP) * T0

KAPPA = 20.0
LEAK_DAMP_PS = 200.0
R_LEAK_STAR = 2.0 * (T0 * 1e12) / LEAK_DAMP_PS


def epstein_gamma(temp_k: float, delta: float = 1.0) -> float:
    cbar = math.sqrt(8.0 * KB * temp_k / (math.pi * M_AR))
    return delta * (4.0 * math.pi / 3.0) * R_NP**2 * RHO_AR * cbar


def interpolate_s(gamma_star: float) -> tuple[list[list[float]], list[list[float]]]:
    """Compute S and Cvv = I-S for a requested dimensionless friction.

    Reuses the covariance routines from oscillatory_reservoir_python by
    temporarily changing its global GAMMA.
    """
    import oscillatory_reservoir_python as model

    old_gamma = model.GAMMA
    model.GAMMA = gamma_star
    try:
        s, c_vv = heat_transfer_matrix()
    finally:
        model.GAMMA = old_gamma
    return s, c_vv


def matvec(a: list[list[float]], x: list[float]) -> list[float]:
    return [sum(a[i][j] * x[j] for j in range(len(x))) for i in range(len(a))]


def rhs(temp: list[float], s: list[list[float]], tamb: float, alpha_heat: float) -> list[float]:
    st = matvec(s, temp)
    return [alpha_heat * KAPPA * GAMMA_STAR * st[i] - R_LEAK_STAR * (temp[i] - tamb) for i in range(N)]


def rk4_step(temp: list[float], dtau: float, s: list[list[float]], tamb: float, alpha_heat: float) -> list[float]:
    k1 = rhs(temp, s, tamb, alpha_heat)
    k2 = rhs([temp[i] + 0.5 * dtau * k1[i] for i in range(N)], s, tamb, alpha_heat)
    k3 = rhs([temp[i] + 0.5 * dtau * k2[i] for i in range(N)], s, tamb, alpha_heat)
    k4 = rhs([temp[i] + dtau * k3[i] for i in range(N)], s, tamb, alpha_heat)
    return [temp[i] + dtau * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0 for i in range(N)]


def simulate_case(
    name: str,
    temps0: list[float],
    tamb: float,
    tmax_ps: float,
    dt_ps: float,
    s: list[list[float]],
    c_vv: list[list[float]],
    alpha_heat: float,
) -> list[list[float]]:
    rows = []
    temp = temps0[:]
    dtau = dt_ps / (T0 * 1e12)
    steps = int(round(tmax_ps / dt_ps))
    for step in range(steps + 1):
        t_ps = step * dt_ps
        teff = matvec(c_vv, temp)
        rows.append([t_ps] + temp[:] + teff[:])
        if step < steps:
            temp = rk4_step(temp, dtau, s, tamb, alpha_heat)

    write_csv(
        OUT / f"{name}_analytic.csv",
        ["time_ps", "Tar1", "Tar2", "Tar3", "Tnp1_analytic", "Tnp2_analytic", "Tnp3_analytic"],
        rows,
    )
    write_svg(
        OUT / f"{name}_reservoirs.svg",
        rows,
        [1, 2, 3],
        ["Tar1 analytic", "Tar2 analytic", "Tar3 analytic"],
        f"{name}: analytic reservoir temperatures",
        "temperature (K)",
    )
    write_svg(
        OUT / f"{name}_nanoparticles.svg",
        rows,
        [4, 5, 6],
        ["Tnp1 analytic", "Tnp2 analytic", "Tnp3 analytic"],
        f"{name}: analytic Cu kinetic temperatures",
        "temperature (K)",
    )
    return rows


def simulate_alpha_scan(
    name: str,
    temps0: list[float],
    tamb: float,
    tmax_ps: float,
    dt_ps: float,
    s: list[list[float]],
    alpha_values: list[float],
) -> None:
    all_rows = []
    for alpha in alpha_values:
        temp = temps0[:]
        dtau = dt_ps / (T0 * 1e12)
        steps = int(round(tmax_ps / dt_ps))
        for step in range(steps + 1):
            if step % 5 == 0:
                all_rows.append([alpha, step * dt_ps] + temp[:])
            if step < steps:
                temp = rk4_step(temp, dtau, s, tamb, alpha)
    write_csv(OUT / f"{name}_alpha_scan.csv", ["alphaHeat", "time_ps", "Tar1", "Tar2", "Tar3"], all_rows)

    # Plot mean reservoir temperature for each alpha. This is the most direct
    # diagnostic for matching the LAMMPS relaxation envelope.
    rows_by_alpha = []
    for alpha in alpha_values:
        curve = [[r[1], sum(r[2:5]) / 3.0] for r in all_rows if abs(r[0] - alpha) < 1e-15]
        rows_by_alpha.append((alpha, curve))
    write_multicurve_svg(
        OUT / f"{name}_alpha_scan_mean_reservoir.svg",
        rows_by_alpha,
        f"{name}: mean reservoir temperature for heat-rate scale alphaHeat",
        "mean reservoir temperature (K)",
    )


def write_csv(path: Path, header: list[str], rows: list[list[float | str]]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def write_svg(path: Path, rows: list[list[float]], cols: list[int], labels: list[str], title: str, ylabel: str) -> None:
    width, height = 1050, 650
    left, right, top, bottom = 85, 190, 60, 70
    xs = [row[0] for row in rows]
    ys = [row[c] for row in rows for c in cols]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    pad = 0.08 * max(1.0, ymax - ymin)
    ymin -= pad
    ymax += pad
    colors = ["#2ca02c", "#1f77b4", "#ff00cc", "#d62728", "#9467bd", "#ff7f0e"]

    def sx(x: float) -> float:
        return left + (x - xmin) / (xmax - xmin) * (width - left - right)

    def sy(y: float) -> float:
        return height - bottom - (y - ymin) / (ymax - ymin) * (height - top - bottom)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width/2}" y="33" text-anchor="middle" font-family="Arial" font-size="20">{title}</text>',
        f'<line x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}" stroke="black"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" stroke="black"/>',
        f'<text x="{width/2}" y="{height-22}" text-anchor="middle" font-family="Arial" font-size="15">time (ps)</text>',
        f'<text x="24" y="{height/2}" transform="rotate(-90 24 {height/2})" text-anchor="middle" font-family="Arial" font-size="15">{ylabel}</text>',
    ]

    for frac in [0.0, 0.25, 0.5, 0.75, 1.0]:
        x = xmin + frac * (xmax - xmin)
        xp = sx(x)
        parts.append(f'<line x1="{xp:.2f}" y1="{top}" x2="{xp:.2f}" y2="{height-bottom}" stroke="#ddd"/>')
        parts.append(f'<text x="{xp:.2f}" y="{height-bottom+22}" text-anchor="middle" font-family="Arial" font-size="12">{x:.0f}</text>')
        y = ymin + frac * (ymax - ymin)
        yp = sy(y)
        parts.append(f'<line x1="{left}" y1="{yp:.2f}" x2="{width-right}" y2="{yp:.2f}" stroke="#ddd"/>')
        parts.append(f'<text x="{left-10}" y="{yp+4:.2f}" text-anchor="end" font-family="Arial" font-size="12">{y:.1f}</text>')

    for idx, c in enumerate(cols):
        pts = " ".join(f"{sx(row[0]):.2f},{sy(row[c]):.2f}" for row in rows)
        color = colors[idx % len(colors)]
        parts.append(f'<polyline points="{pts}" fill="none" stroke="{color}" stroke-width="2.2"/>')
        yleg = top + 25 * idx
        parts.append(f'<line x1="{width-right+20}" y1="{yleg}" x2="{width-right+55}" y2="{yleg}" stroke="{color}" stroke-width="3"/>')
        parts.append(f'<text x="{width-right+62}" y="{yleg+5}" font-family="Arial" font-size="13">{labels[idx]}</text>')

    parts.append("</svg>")
    path.write_text("\n".join(parts))


def write_multicurve_svg(path: Path, curves: list[tuple[float, list[list[float]]]], title: str, ylabel: str) -> None:
    width, height = 1050, 650
    left, right, top, bottom = 85, 220, 60, 70
    xs = [p[0] for _, curve in curves for p in curve]
    ys = [p[1] for _, curve in curves for p in curve]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    pad = 0.08 * max(1.0, ymax - ymin)
    ymin -= pad
    ymax += pad
    colors = ["#d62728", "#ff7f0e", "#2ca02c", "#1f77b4", "#9467bd"]

    def sx(x: float) -> float:
        return left + (x - xmin) / (xmax - xmin) * (width - left - right)

    def sy(y: float) -> float:
        return height - bottom - (y - ymin) / (ymax - ymin) * (height - top - bottom)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width/2}" y="33" text-anchor="middle" font-family="Arial" font-size="20">{title}</text>',
        f'<line x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}" stroke="black"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" stroke="black"/>',
        f'<text x="{width/2}" y="{height-22}" text-anchor="middle" font-family="Arial" font-size="15">time (ps)</text>',
        f'<text x="24" y="{height/2}" transform="rotate(-90 24 {height/2})" text-anchor="middle" font-family="Arial" font-size="15">{ylabel}</text>',
    ]
    for frac in [0.0, 0.25, 0.5, 0.75, 1.0]:
        x = xmin + frac * (xmax - xmin)
        xp = sx(x)
        parts.append(f'<line x1="{xp:.2f}" y1="{top}" x2="{xp:.2f}" y2="{height-bottom}" stroke="#ddd"/>')
        parts.append(f'<text x="{xp:.2f}" y="{height-bottom+22}" text-anchor="middle" font-family="Arial" font-size="12">{x:.0f}</text>')
        y = ymin + frac * (ymax - ymin)
        yp = sy(y)
        parts.append(f'<line x1="{left}" y1="{yp:.2f}" x2="{width-right}" y2="{yp:.2f}" stroke="#ddd"/>')
        parts.append(f'<text x="{left-10}" y="{yp+4:.2f}" text-anchor="end" font-family="Arial" font-size="12">{y:.1f}</text>')

    for idx, (alpha, curve) in enumerate(curves):
        pts = " ".join(f"{sx(x):.2f},{sy(y):.2f}" for x, y in curve)
        color = colors[idx % len(colors)]
        parts.append(f'<polyline points="{pts}" fill="none" stroke="{color}" stroke-width="2.2"/>')
        yleg = top + 25 * idx
        parts.append(f'<line x1="{width-right+20}" y1="{yleg}" x2="{width-right+55}" y2="{yleg}" stroke="{color}" stroke-width="3"/>')
        parts.append(f'<text x="{width-right+62}" y="{yleg+5}" font-family="Arial" font-size="13">alpha={alpha:g}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts))


def make_tables(s: list[list[float]], c_vv: list[list[float]]) -> str:
    ep_rows = []
    for temp in [110.0, 140.0, 170.0]:
        g1 = epstein_gamma(temp, 1.0)
        g139 = epstein_gamma(temp, 1.39)
        ep_rows.append([temp, g1, M_NP / g1 * 1e12, g139, M_NP / g139 * 1e12])

    write_csv(
        OUT / "epstein_gamma.csv",
        ["T_K", "Gamma_delta1_kg_s", "M_over_Gamma_delta1_ps", "Gamma_delta1p39_kg_s", "M_over_Gamma_delta1p39_ps"],
        ep_rows,
    )

    s_eigs = eig3(s)
    a_res = [[KAPPA * GAMMA_STAR * s[i][j] - (R_LEAK_STAR if i == j else 0.0) for j in range(N)] for i in range(N)]
    a_eigs = eig3(a_res)

    summary_rows = [
        ["N_Cu per nanoparticle", N_CU, "atoms"],
        ["M_NP", M_NP, "kg"],
        ["R_NP", R_NP * 1e9, "nm"],
        ["Ar number density", N_DENSITY_AR_NM3, "nm^-3"],
        ["Ar mass density", RHO_AR, "kg/m^3"],
        ["K diagonal", KFORCE_EV_A2, "eV/A^2"],
        ["K diagonal", K_DIAG, "N/m"],
        ["t0=sqrt(M/K)", T0 * 1e12, "ps"],
        ["Gamma_ref", GAMMA_REF, "kg/s"],
        ["Gamma_star", GAMMA_STAR, "dimensionless"],
        ["M/Gamma_ref", M_NP / GAMMA_REF * 1e12, "ps"],
        ["M/(2 Gamma_ref)", M_NP / (2.0 * GAMMA_REF) * 1e12, "ps"],
        ["leakDamp", LEAK_DAMP_PS, "ps"],
        ["r_star=2 t0/leakDamp", R_LEAK_STAR, "dimensionless"],
    ]
    write_csv(OUT / "mapping_summary.csv", ["quantity", "value", "unit"], summary_rows)

    md = []
    md.append("# Analytic LAMMPS Mapping Tables\n")
    md.append("## Physical Parameters\n")
    md.append("| quantity | value | unit |")
    md.append("|---|---:|---|")
    for q, v, u in summary_rows:
        md.append(f"| {q} | {float(v):.6g} | {u} |")

    md.append("\n## Epstein Drag Estimates\n")
    md.append("| T (K) | Gamma delta=1 (kg/s) | M/Gamma (ps) | Gamma delta=1.39 (kg/s) | M/Gamma (ps) |")
    md.append("|---:|---:|---:|---:|---:|")
    for row in ep_rows:
        md.append(f"| {row[0]:.0f} | {row[1]:.4e} | {row[2]:.2f} | {row[3]:.4e} | {row[4]:.2f} |")

    md.append("\n## Heat Transfer Matrix S for Current Physical Mapping\n")
    md.append("Rows define `J = gamma S T` in the dimensionless oscillator model.")
    md.append("")
    md.append("| row | col1 | col2 | col3 |")
    md.append("|---:|---:|---:|---:|")
    for i, row in enumerate(s, start=1):
        md.append(f"| {i} | {row[0]:.6g} | {row[1]:.6g} | {row[2]:.6g} |")

    md.append("\n## Effective Cu Kinetic Temperature Response Cvv = I - S\n")
    md.append("Analytic nanoparticle kinetic temperatures are `Tnp = Cvv . Tar`.")
    md.append("")
    md.append("| row | col1 | col2 | col3 |")
    md.append("|---:|---:|---:|---:|")
    for i, row in enumerate(c_vv, start=1):
        md.append(f"| {i} | {row[0]:.6g} | {row[1]:.6g} | {row[2]:.6g} |")

    md.append("\n## Eigenvalues\n")
    md.append("| matrix | eigenvalue |")
    md.append("|---|---:|")
    for z in s_eigs:
        md.append(f"| S | {z.real:.6g} {z.imag:+.6g} i |")
    for z in a_eigs:
        period = "" if abs(z.imag) < 1e-12 else f", period={2.0*math.pi/abs(z.imag)*T0*1e12:.1f} ps"
        decay = "" if z.real >= 0 else f", decay={T0*1e12/abs(z.real):.1f} ps"
        md.append(f"| reservoir A | {z.real:.6g} {z.imag:+.6g} i{period}{decay} |")

    md.append("\n## Heat-Rate Scale alphaHeat\n")
    md.append("The toy analytic value `kappa=20` is kept, but the atomistic reservoir has a large finite heat capacity.")
    md.append("Use `alphaHeat` as a calibration multiplier in")
    md.append("")
    md.append("`dT/dtau = alphaHeat*kappa*GammaStar*S*T - rStar*(T-Tamb)`.")
    md.append("")
    md.append("| alphaHeat | final mean T at 1800 ps for current LAMMPS initial condition (K) |")
    md.append("|---:|---:|")
    for alpha in [1.0, 0.3, 0.1, 0.03, 0.01]:
        temp = [173.5, 143.5, 125.5]
        dtau = 1.0 / (T0 * 1e12)
        for _ in range(1800):
            temp = rk4_step(temp, dtau, s, 110.0, alpha)
        md.append(f"| {alpha:g} | {sum(temp)/3.0:.3f} |")

    text = "\n".join(md) + "\n"
    (OUT / "parameter_tables.md").write_text(text)
    return text


def main() -> None:
    OUT.mkdir(exist_ok=True)
    s, c_vv = interpolate_s(GAMMA_STAR)
    tables = make_tables(s, c_vv)
    simulate_case("current_lammps", [173.5, 143.5, 125.5], 110.0, 1800.0, 1.0, s, c_vv, 1.0)
    simulate_case("centered_excitation", [170.0, 140.0, 110.0], 140.0, 3000.0, 1.0, s, c_vv, 1.0)
    simulate_alpha_scan("current_lammps", [173.5, 143.5, 125.5], 110.0, 1800.0, 1.0, s, [1.0, 0.3, 0.1, 0.03, 0.01])

    print(tables)
    print(f"Wrote analytic results to {OUT.resolve()}")


if __name__ == "__main__":
    main()
