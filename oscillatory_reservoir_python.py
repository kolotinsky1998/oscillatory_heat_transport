#!/usr/bin/env python3
"""Python demonstration of oscillatory heat transport through 3 oscillators.

No third-party packages are required.  The script:

1. Builds the heat-transfer matrix S for the force matrix F used in the
   Mathematica notebook/script.
2. Checks the complex spectrum of S and the stabilized reservoir matrix.
3. Compares the reservoir-temperature ODE against the closed analytic curve.
4. Runs an underdamped Langevin ensemble for the same nonreciprocal oscillator
   system and compares stationary kinetic temperatures with the covariance
   analytics.
5. Writes CSV and SVG plots.
"""

from __future__ import annotations

import cmath
import csv
import math
import random
from pathlib import Path


N = 3
GAMMA = 0.1
KAPPA = 20.0
DRAIN = 0.98
TAMB = 1.0

F = [
    [1.00, 0.14, 0.06],
    [0.06, 1.00, 0.14],
    [0.14, 0.06, 1.00],
]

T0 = [1.45, 0.95, 0.65]

OUTDIR = Path("python_results")


def solve_linear(a: list[list[complex]], b: list[complex]) -> list[complex]:
    """Dense Gaussian elimination with partial pivoting."""
    n = len(b)
    m = [list(map(complex, a[i])) + [complex(b[i])] for i in range(n)]

    for k in range(n):
        pivot = max(range(k, n), key=lambda r: abs(m[r][k]))
        if abs(m[pivot][k]) < 1e-14:
            raise ValueError("singular linear system")
        m[k], m[pivot] = m[pivot], m[k]

        pv = m[k][k]
        for j in range(k, n + 1):
            m[k][j] /= pv

        for i in range(n):
            if i == k:
                continue
            factor = m[i][k]
            if factor == 0:
                continue
            for j in range(k, n + 1):
                m[i][j] -= factor * m[k][j]

    return [m[i][n] for i in range(n)]


def matmul(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [
        [sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0]))]
        for i in range(len(a))
    ]


def matvec(a: list[list[float]], x: list[float]) -> list[float]:
    return [sum(a[i][j] * x[j] for j in range(len(x))) for i in range(len(a))]


def add_vec(a: list[float], b: list[float]) -> list[float]:
    return [a[i] + b[i] for i in range(len(a))]


def scale_vec(s: float, a: list[float]) -> list[float]:
    return [s * x for x in a]


def eig3(m: list[list[float]]) -> list[complex]:
    """Eigenvalues of a real 3x3 matrix via Durand-Kerner cubic roots."""
    tr = sum(m[i][i] for i in range(3))
    m2 = matmul(m, m)
    tr2 = sum(m2[i][i] for i in range(3))
    det = (
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )

    # det(lambda I - M) = lambda^3 + a lambda^2 + b lambda + c
    a = -tr
    b = 0.5 * (tr * tr - tr2)
    c = -det

    roots = [1 + 0j, -0.5 + 0.8j, -0.5 - 0.8j]

    def p(z: complex) -> complex:
        return z**3 + a * z**2 + b * z + c

    for _ in range(200):
        new_roots = []
        for i, z in enumerate(roots):
            denom = 1 + 0j
            for j, w in enumerate(roots):
                if i != j:
                    denom *= z - w
            new_roots.append(z - p(z) / denom)
        if max(abs(new_roots[i] - roots[i]) for i in range(3)) < 1e-13:
            roots = new_roots
            break
        roots = new_roots
    return roots


def covariance_for_temperatures(f: list[list[float]], gamma: float, temps: list[float]) -> list[list[float]]:
    """Solve A C + C A^T + Q = 0 for x=(q1,q2,q3,v1,v2,v3)."""
    n = len(f)
    dim = 2 * n
    a = [[0.0 for _ in range(dim)] for _ in range(dim)]
    for i in range(n):
        a[i][n + i] = 1.0
        for j in range(n):
            a[n + i][j] = -f[i][j]
        a[n + i][n + i] = -gamma

    system = []
    rhs = []
    for i in range(dim):
        for j in range(dim):
            row = []
            for alpha in range(dim):
                for beta in range(dim):
                    value = 0.0
                    if beta == j:
                        value += a[i][alpha]
                    if alpha == i:
                        value += a[j][beta]
                    row.append(value)
            qij = 0.0
            if i >= n and j >= n and i == j:
                qij = 2.0 * gamma * temps[i - n]
            system.append(row)
            rhs.append(-qij)

    solution = solve_linear(system, rhs)
    cov = [[0.0 for _ in range(dim)] for _ in range(dim)]
    for i in range(dim):
        for j in range(dim):
            cov[i][j] = solution[i * dim + j].real
    return cov


def effective_temperature_response() -> list[list[float]]:
    response = [[0.0 for _ in range(N)] for _ in range(N)]
    for j in range(N):
        temps = [0.0] * N
        temps[j] = 1.0
        cov = covariance_for_temperatures(F, GAMMA, temps)
        for i in range(N):
            response[i][j] = cov[N + i][N + i]
    return response


def heat_transfer_matrix() -> tuple[list[list[float]], list[list[float]]]:
    c_vv = effective_temperature_response()
    s = [[(1.0 if i == j else 0.0) - c_vv[i][j] for j in range(N)] for i in range(N)]
    return s, c_vv


def reservoir_matrix(s: list[list[float]]) -> list[list[float]]:
    return [
        [KAPPA * GAMMA * s[i][j] - (DRAIN if i == j else 0.0) for j in range(N)]
        for i in range(N)
    ]


def reservoir_rhs(a_res: list[list[float]], t: list[float]) -> list[float]:
    return add_vec(matvec(a_res, t), [DRAIN * TAMB] * N)


def rk4_step(a_res: list[list[float]], temp: list[float], dt: float) -> list[float]:
    k1 = reservoir_rhs(a_res, temp)
    k2 = reservoir_rhs(a_res, add_vec(temp, scale_vec(0.5 * dt, k1)))
    k3 = reservoir_rhs(a_res, add_vec(temp, scale_vec(0.5 * dt, k2)))
    k4 = reservoir_rhs(a_res, add_vec(temp, scale_vec(dt, k3)))
    return [
        temp[i] + dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0
        for i in range(N)
    ]


def analytic_reservoir_temperature(time: float) -> list[float]:
    """Closed-form solution for the circulant matrix used here."""
    # These constants are calculated from A = kappa gamma S - drain I.
    tss = 0.39140357578686447
    cfast = 0.6252630908798021
    lamfast = -2.503809522
    lamr = -0.042796722
    omega = 0.13966544107302703
    aa = [0.43333333333333324, -0.06666666666666665, -0.3666666666666666]
    bb = [-0.17320508075688648, 0.461880215351701, -0.28867513459481287]
    return [
        tss
        + cfast * math.exp(lamfast * time)
        + math.exp(lamr * time) * (aa[i] * math.cos(omega * time) + bb[i] * math.sin(omega * time))
        for i in range(N)
    ]


def simulate_stationary_langevin(
    temps: list[float],
    dt: float = 0.005,
    steps: int = 30_000,
    burn: int = 5_000,
    ntraj: int = 250,
    sample_stride: int = 10,
    seed: int = 91,
) -> list[float]:
    """Ensemble Langevin check for constant reservoirs.

    With gamma=0.1 the correlation time is long enough that a single trajectory
    gives poor v^2 statistics.  Many modest trajectories are faster and much
    less misleading in pure Python.
    """
    rng = random.Random(seed)
    q = [[0.1 * rng.gauss(0.0, 1.0) for _ in range(N)] for _ in range(ntraj)]
    v = [[math.sqrt(max(0.0, temps[i])) * rng.gauss(0.0, 1.0) for i in range(N)] for _ in range(ntraj)]
    ou = math.exp(-GAMMA * dt)
    noise = [math.sqrt(max(0.0, temps[i]) * (1.0 - ou * ou)) for i in range(N)]
    sums = [0.0, 0.0, 0.0]
    count = 0

    def force(pos: list[float]) -> list[float]:
        return [-sum(F[i][j] * pos[j] for j in range(N)) for i in range(N)]

    for step in range(steps):
        for r in range(ntraj):
            frc = force(q[r])
            for i in range(N):
                v[r][i] += 0.5 * dt * frc[i]
                v[r][i] = ou * v[r][i] + noise[i] * rng.gauss(0.0, 1.0)
                q[r][i] += dt * v[r][i]
            frc = force(q[r])
            for i in range(N):
                v[r][i] += 0.5 * dt * frc[i]

        if step >= burn and step % sample_stride == 0:
            for r in range(ntraj):
                for i in range(N):
                    sums[i] += v[r][i] * v[r][i]
            count += ntraj

    return [s / count for s in sums]


def simulate_dynamic_langevin(
    dt: float = 0.01, tmax: float = 160.0, ntraj: int = 400, sample_every: int = 20, seed: int = 123
) -> list[list[float]]:
    """Ensemble Langevin run driven by the analytic reservoir temperatures."""
    rng = random.Random(seed)
    steps = int(round(tmax / dt))
    q = [[0.0, 0.0, 0.0] for _ in range(ntraj)]
    v = [[0.0, 0.0, 0.0] for _ in range(ntraj)]
    for r in range(ntraj):
        q[r] = [0.35 + 0.01 * rng.gauss(0.0, 1.0), -0.10, -0.30]

    ou = math.exp(-GAMMA * dt)
    rows = []

    def forces(pos: list[float]) -> list[float]:
        return [-sum(F[i][j] * pos[j] for j in range(N)) for i in range(N)]

    for step in range(steps + 1):
        time = step * dt
        target = analytic_reservoir_temperature(time)
        if step % sample_every == 0:
            kin = [0.0, 0.0, 0.0]
            mean_q = [0.0, 0.0, 0.0]
            for r in range(ntraj):
                for i in range(N):
                    kin[i] += v[r][i] * v[r][i]
                    mean_q[i] += q[r][i]
            rows.append(
                [time]
                + target
                + [kin[i] / ntraj for i in range(N)]
                + [mean_q[i] / ntraj for i in range(N)]
            )
        if step == steps:
            break

        sigma = [math.sqrt(max(0.0, target[i]) * (1.0 - ou * ou)) for i in range(N)]
        for r in range(ntraj):
            frc = forces(q[r])
            for i in range(N):
                v[r][i] += 0.5 * dt * frc[i]
                v[r][i] = ou * v[r][i] + sigma[i] * rng.gauss(0.0, 1.0)
                q[r][i] += dt * v[r][i]
            frc = forces(q[r])
            for i in range(N):
                v[r][i] += 0.5 * dt * frc[i]

    return rows


def write_csv(path: Path, header: list[str], rows: list[list[float]]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def write_svg(path: Path, rows: list[list[float]], cols: list[int], labels: list[str], title: str, ylabel: str) -> None:
    width, height = 1000, 620
    left, right, top, bottom = 80, 170, 55, 70
    xs = [row[0] for row in rows]
    ys = [row[c] for row in rows for c in cols]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    pad = 0.08 * (ymax - ymin if ymax > ymin else 1.0)
    ymin -= pad
    ymax += pad
    colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#17becf"]

    def sx(x: float) -> float:
        return left + (x - xmin) / (xmax - xmin) * (width - left - right)

    def sy(y: float) -> float:
        return height - bottom - (y - ymin) / (ymax - ymin) * (height - top - bottom)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width/2}" y="30" text-anchor="middle" font-family="Arial" font-size="20">{title}</text>',
        f'<line x1="{left}" y1="{height-bottom}" x2="{width-right}" y2="{height-bottom}" stroke="black"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{height-bottom}" stroke="black"/>',
        f'<text x="{width/2}" y="{height-20}" text-anchor="middle" font-family="Arial" font-size="15">time</text>',
        f'<text x="22" y="{height/2}" transform="rotate(-90 22 {height/2})" text-anchor="middle" font-family="Arial" font-size="15">{ylabel}</text>',
    ]

    for frac in [0, 0.25, 0.5, 0.75, 1.0]:
        x = xmin + frac * (xmax - xmin)
        xp = sx(x)
        parts.append(f'<line x1="{xp:.2f}" y1="{top}" x2="{xp:.2f}" y2="{height-bottom}" stroke="#ddd"/>')
        parts.append(f'<text x="{xp:.2f}" y="{height-bottom+22}" text-anchor="middle" font-family="Arial" font-size="12">{x:.0f}</text>')
        y = ymin + frac * (ymax - ymin)
        yp = sy(y)
        parts.append(f'<line x1="{left}" y1="{yp:.2f}" x2="{width-right}" y2="{yp:.2f}" stroke="#ddd"/>')
        parts.append(f'<text x="{left-10}" y="{yp+4:.2f}" text-anchor="end" font-family="Arial" font-size="12">{y:.3g}</text>')

    for idx, c in enumerate(cols):
        pts = " ".join(f"{sx(row[0]):.2f},{sy(row[c]):.2f}" for row in rows)
        color = colors[idx % len(colors)]
        parts.append(f'<polyline points="{pts}" fill="none" stroke="{color}" stroke-width="2"/>')
        yleg = top + 24 * idx
        parts.append(f'<line x1="{width-right+25}" y1="{yleg}" x2="{width-right+55}" y2="{yleg}" stroke="{color}" stroke-width="3"/>')
        parts.append(f'<text x="{width-right+62}" y="{yleg+5}" font-family="Arial" font-size="13">{labels[idx]}</text>')

    parts.append("</svg>")
    path.write_text("\n".join(parts))


def main() -> None:
    OUTDIR.mkdir(exist_ok=True)

    s, c_vv = heat_transfer_matrix()
    s_eigs = eig3(s)
    a_res = reservoir_matrix(s)
    a_eigs = eig3(a_res)

    print("F matrix:")
    for row in F:
        print("  ", " ".join(f"{x: .8f}" for x in row))
    print("\nS matrix, J = gamma S T:")
    for row in s:
        print("  ", " ".join(f"{x: .8f}" for x in row))
    print("\neig(S):")
    for z in s_eigs:
        print(f"  {z.real: .10f} {z.imag:+.10f}i")
    print("\neig(kappa gamma S - drain I):")
    for z in a_eigs:
        print(f"  {z.real: .10f} {z.imag:+.10f}i")

    # Reservoir ODE: RK4 versus analytic closed form.
    dt = 0.02
    tmax = 160.0
    steps = int(round(tmax / dt))
    temp = T0[:]
    reservoir_rows = []
    max_err = 0.0
    for step in range(steps + 1):
        time = step * dt
        exact = analytic_reservoir_temperature(time)
        max_err = max(max_err, max(abs(temp[i] - exact[i]) for i in range(N)))
        if step % 10 == 0:
            reservoir_rows.append([time] + temp[:] + exact[:] + [temp[i] - exact[i] for i in range(N)])
        if step < steps:
            temp = rk4_step(a_res, temp, dt)

    print(f"\nRK4 reservoir max error against closed analytic curve: {max_err:.3e}")

    # Stationary Langevin check against covariance analytics.
    fixed_t = T0
    cov = covariance_for_temperatures(F, GAMMA, fixed_t)
    analytic_teff = [cov[N + i][N + i] for i in range(N)]
    measured_teff = simulate_stationary_langevin(fixed_t)
    print("\nStationary Langevin check for fixed T = [1.45, 0.95, 0.65]:")
    for i in range(N):
        rel = (measured_teff[i] - analytic_teff[i]) / analytic_teff[i]
        print(
            f"  <v{i+1}^2>: simulation={measured_teff[i]:.5f}, "
            f"analytic={analytic_teff[i]:.5f}, rel.err={rel:+.2%}"
        )

    # Dynamic ensemble driven by the oscillatory reservoir temperatures.
    dynamic_rows = simulate_dynamic_langevin()

    write_csv(
        OUTDIR / "reservoir_rk4_vs_analytic.csv",
        ["time", "T1_rk4", "T2_rk4", "T3_rk4", "T1_exact", "T2_exact", "T3_exact", "dT1", "dT2", "dT3"],
        reservoir_rows,
    )
    write_csv(
        OUTDIR / "dynamic_langevin.csv",
        ["time", "T1", "T2", "T3", "kin1", "kin2", "kin3", "mean_x1", "mean_x2", "mean_x3"],
        dynamic_rows,
    )
    write_svg(
        OUTDIR / "reservoir_temperatures.svg",
        reservoir_rows,
        [1, 2, 3],
        ["T1 RK4", "T2 RK4", "T3 RK4"],
        "Reservoir temperatures: damped oscillations",
        "temperature",
    )
    write_svg(
        OUTDIR / "rk4_minus_analytic.svg",
        reservoir_rows,
        [7, 8, 9],
        ["T1 error", "T2 error", "T3 error"],
        "Numerical reservoir ODE minus analytic curve",
        "error",
    )
    write_svg(
        OUTDIR / "dynamic_langevin_kinetic.svg",
        dynamic_rows,
        [1, 2, 3, 4, 5, 6],
        ["bath T1", "bath T2", "bath T3", "sim <v1^2>", "sim <v2^2>", "sim <v3^2>"],
        "Driven Langevin ensemble: bath targets and kinetic temperatures",
        "temperature",
    )

    print(f"\nWrote results to: {OUTDIR.resolve()}")
    print("  reservoir_rk4_vs_analytic.csv")
    print("  dynamic_langevin.csv")
    print("  reservoir_temperatures.svg")
    print("  rk4_minus_analytic.svg")
    print("  dynamic_langevin_kinetic.svg")


if __name__ == "__main__":
    main()
