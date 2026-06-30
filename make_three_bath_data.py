#!/usr/bin/env python3
"""Generate a three-reservoir Cu/Ar LAMMPS data file.

The geometry follows the all-atom spirit of Loos et al., Sci. Rep. 13, 4517
(2023): independent argon reservoirs, copper nanoparticles, and fixed copper
walls.  The resulting data file is used by in.three_bath_oscillatory_heat.lmp.

No third-party Python packages are required.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path


def fcc_basis(a: float) -> list[tuple[float, float, float]]:
    return [
        (0.0, 0.0, 0.0),
        (0.5 * a, 0.5 * a, 0.0),
        (0.5 * a, 0.0, 0.5 * a),
        (0.0, 0.5 * a, 0.5 * a),
    ]


def add_fcc_box(
    atoms: list[tuple[int, float, float, float]],
    atom_type: int,
    a: float,
    xlo: float,
    xhi: float,
    ylo: float,
    yhi: float,
    zlo: float,
    zhi: float,
    exclude_spheres: list[tuple[float, float, float, float]] | None = None,
) -> None:
    exclude_spheres = exclude_spheres or []
    nx = int(math.ceil((xhi - xlo) / a)) + 1
    ny = int(math.ceil((yhi - ylo) / a)) + 1
    nz = int(math.ceil((zhi - zlo) / a)) + 1
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                base = (xlo + ix * a, ylo + iy * a, zlo + iz * a)
                for bx, by, bz in fcc_basis(a):
                    x = base[0] + bx
                    y = base[1] + by
                    z = base[2] + bz
                    if not (xlo <= x < xhi and ylo <= y < yhi and zlo <= z < zhi):
                        continue
                    excluded = False
                    for cx, cy, cz, radius in exclude_spheres:
                        if (x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2 < radius**2:
                            excluded = True
                            break
                    if not excluded:
                        atoms.append((atom_type, x, y, z))


def add_fcc_sphere(
    atoms: list[tuple[int, float, float, float]],
    atom_type: int,
    a: float,
    cx: float,
    cy: float,
    cz: float,
    radius: float,
) -> int:
    before = len(atoms)
    n = int(math.ceil(radius / a)) + 2
    for ix in range(-n, n + 1):
        for iy in range(-n, n + 1):
            for iz in range(-n, n + 1):
                base = (cx + ix * a, cy + iy * a, cz + iz * a)
                for bx, by, bz in fcc_basis(a):
                    x = base[0] + bx
                    y = base[1] + by
                    z = base[2] + bz
                    if (x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2 <= radius**2:
                        atoms.append((atom_type, x, y, z))
    return len(atoms) - before


def add_wall_slab(
    atoms: list[tuple[int, float, float, float]],
    atom_type: int,
    a: float,
    xcenter: float,
    thickness: float,
    ylo: float,
    yhi: float,
    zlo: float,
    zhi: float,
) -> None:
    add_fcc_box(
        atoms,
        atom_type,
        a,
        xcenter - 0.5 * thickness,
        xcenter + 0.5 * thickness,
        ylo,
        yhi,
        zlo,
        zhi,
    )


def write_data(path: Path, atoms: list[tuple[int, float, float, float]], box: tuple[float, float, float, float, float, float]) -> None:
    xlo, xhi, ylo, yhi, zlo, zhi = box
    with path.open("w") as f:
        f.write("Three independent Ar reservoirs with Cu nanoparticles\n\n")
        f.write(f"{len(atoms)} atoms\n")
        f.write("7 atom types\n\n")
        f.write(f"{xlo:.8f} {xhi:.8f} xlo xhi\n")
        f.write(f"{ylo:.8f} {yhi:.8f} ylo yhi\n")
        f.write(f"{zlo:.8f} {zhi:.8f} zlo zhi\n\n")
        f.write("Masses\n\n")
        f.write("1 39.948  # Ar reservoir 1\n")
        f.write("2 39.948  # Ar reservoir 2\n")
        f.write("3 39.948  # Ar reservoir 3\n")
        f.write("4 63.546  # Cu nanoparticle 1\n")
        f.write("5 63.546  # Cu nanoparticle 2\n")
        f.write("6 63.546  # Cu nanoparticle 3\n")
        f.write("7 63.546  # fixed Cu walls\n\n")
        f.write("Atoms # atomic\n\n")
        for idx, (atom_type, x, y, z) in enumerate(atoms, start=1):
            f.write(f"{idx} {atom_type} {x:.8f} {y:.8f} {z:.8f}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="three_bath_cu_ar.data")
    parser.add_argument("--compartment-length", type=float, default=95.0, help="A")
    parser.add_argument("--box-y", type=float, default=95.0, help="A")
    parser.add_argument("--box-z", type=float, default=95.0, help="A")
    parser.add_argument("--ar-lattice", type=float, default=5.30, help="A, liquid-like initial spacing")
    parser.add_argument("--cu-lattice", type=float, default=3.615, help="A")
    parser.add_argument("--np-radius", type=float, default=14.4, help="A")
    parser.add_argument("--wall-thickness", type=float, default=7.23, help="A")
    parser.add_argument("--wall-margin", type=float, default=9.0, help="A")
    args = parser.parse_args()

    atoms: list[tuple[int, float, float, float]] = []
    lx = args.compartment_length
    ly = args.box_y
    lz = args.box_z
    xlo, xhi = -0.5 * args.wall_thickness, 3.0 * lx + 0.5 * args.wall_thickness
    ylo, yhi = -0.5 * ly, 0.5 * ly
    zlo, zhi = -0.5 * lz, 0.5 * lz

    centers = [
        (0.5 * lx, 0.0, 0.0),
        (1.5 * lx, 0.0, 0.0),
        (2.5 * lx, 0.0, 0.0),
    ]

    # Fixed Cu walls at reservoir boundaries and between reservoirs.
    for wall_x in [0.0, lx, 2.0 * lx, 3.0 * lx]:
        add_wall_slab(
            atoms,
            7,
            args.cu_lattice,
            wall_x,
            args.wall_thickness,
            ylo,
            yhi,
            zlo,
            zhi,
        )

    # Copper nanoparticles.
    np_counts = []
    for i, (cx, cy, cz) in enumerate(centers):
        np_counts.append(add_fcc_sphere(atoms, 4 + i, args.cu_lattice, cx, cy, cz, args.np_radius))

    # Argon reservoirs.  Ar is kept away from walls and from the nanoparticle.
    for i, (cx, cy, cz) in enumerate(centers):
        cxlo = i * lx + args.wall_margin
        cxhi = (i + 1) * lx - args.wall_margin
        add_fcc_box(
            atoms,
            1 + i,
            args.ar_lattice,
            cxlo,
            cxhi,
            ylo + args.wall_margin,
            yhi - args.wall_margin,
            zlo + args.wall_margin,
            zhi - args.wall_margin,
            exclude_spheres=[(cx, cy, cz, args.np_radius + 4.0)],
        )

    write_data(Path(args.output), atoms, (xlo, xhi, ylo, yhi, zlo, zhi))
    counts = {t: 0 for t in range(1, 8)}
    for atom_type, *_ in atoms:
        counts[atom_type] += 1
    print(f"Wrote {args.output}")
    print(f"Total atoms: {len(atoms)}")
    print(f"Cu nanoparticle atom counts: {np_counts}")
    for atom_type in range(1, 8):
        print(f"type {atom_type}: {counts[atom_type]}")


if __name__ == "__main__":
    main()
