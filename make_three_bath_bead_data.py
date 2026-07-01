#!/usr/bin/env python3
"""Generate a small three-reservoir Ar/bead LAMMPS data file.

This is a deliberately coarse-grained model for observing finite-reservoir
temperature redistribution in direct MD.  Each reservoir contains a small
argon bath and one heavy bead oscillator.  Fixed wall atoms separate the
reservoirs and confine the gas.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path


Atom = tuple[int, float, float, float]


def fcc_basis(a: float) -> list[tuple[float, float, float]]:
    return [
        (0.0, 0.0, 0.0),
        (0.5 * a, 0.5 * a, 0.0),
        (0.5 * a, 0.0, 0.5 * a),
        (0.0, 0.5 * a, 0.5 * a),
    ]


def add_fcc_box(
    atoms: list[Atom],
    atom_type: int,
    a: float,
    xlo: float,
    xhi: float,
    ylo: float,
    yhi: float,
    zlo: float,
    zhi: float,
    limit: int | None = None,
    exclude_spheres: list[tuple[float, float, float, float]] | None = None,
) -> None:
    exclude_spheres = exclude_spheres or []
    nx = int(math.ceil((xhi - xlo) / a)) + 1
    ny = int(math.ceil((yhi - ylo) / a)) + 1
    nz = int(math.ceil((zhi - zlo) / a)) + 1
    added = 0
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
                    if any((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2 < r**2 for cx, cy, cz, r in exclude_spheres):
                        continue
                    atoms.append((atom_type, x, y, z))
                    added += 1
                    if limit is not None and added >= limit:
                        return


def write_data(path: Path, atoms: list[Atom], box: tuple[float, float, float, float, float, float], bead_mass: float) -> None:
    xlo, xhi, ylo, yhi, zlo, zhi = box
    with path.open("w") as f:
        f.write("Coarse-grained three Ar reservoirs with heavy bead oscillators\n\n")
        f.write(f"{len(atoms)} atoms\n")
        f.write("7 atom types\n\n")
        f.write(f"{xlo:.8f} {xhi:.8f} xlo xhi\n")
        f.write(f"{ylo:.8f} {yhi:.8f} ylo yhi\n")
        f.write(f"{zlo:.8f} {zhi:.8f} zlo zhi\n\n")
        f.write("Masses\n\n")
        f.write("1 39.948  # Ar reservoir 1\n")
        f.write("2 39.948  # Ar reservoir 2\n")
        f.write("3 39.948  # Ar reservoir 3\n")
        f.write(f"4 {bead_mass:.8f}  # coarse bead 1\n")
        f.write(f"5 {bead_mass:.8f}  # coarse bead 2\n")
        f.write(f"6 {bead_mass:.8f}  # coarse bead 3\n")
        f.write("7 63.546  # fixed wall atoms\n\n")
        f.write("Atoms # atomic\n\n")
        for idx, (atom_type, x, y, z) in enumerate(atoms, start=1):
            f.write(f"{idx} {atom_type} {x:.8f} {y:.8f} {z:.8f}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="three_bath_bead.data")
    parser.add_argument("--ar-per-reservoir", type=int, default=200)
    parser.add_argument("--bead-mass", type=float, default=8000.0, help="amu")
    parser.add_argument("--compartment-length", type=float, default=42.0, help="A")
    parser.add_argument("--box-y", type=float, default=30.0, help="A")
    parser.add_argument("--box-z", type=float, default=30.0, help="A")
    parser.add_argument("--ar-lattice", type=float, default=5.2, help="A")
    parser.add_argument("--wall-lattice", type=float, default=3.615, help="A")
    parser.add_argument("--wall-thickness", type=float, default=3.615, help="A")
    parser.add_argument("--wall-margin", type=float, default=4.0, help="A")
    args = parser.parse_args()

    atoms: list[Atom] = []
    lx = args.compartment_length
    ly = args.box_y
    lz = args.box_z
    xlo = -0.5 * args.wall_thickness
    xhi = 3.0 * lx + 0.5 * args.wall_thickness
    ylo, yhi = -0.5 * ly, 0.5 * ly
    zlo, zhi = -0.5 * lz, 0.5 * lz

    centers = [
        (0.5 * lx, 0.0, 0.0),
        (1.5 * lx, 0.0, 0.0),
        (2.5 * lx, 0.0, 0.0),
    ]

    # Fixed wall slabs at the left edge, between reservoirs, and at the right edge.
    for wall_x in [0.0, lx, 2.0 * lx, 3.0 * lx]:
        add_fcc_box(
            atoms,
            7,
            args.wall_lattice,
            wall_x - 0.5 * args.wall_thickness,
            wall_x + 0.5 * args.wall_thickness,
            ylo,
            yhi,
            zlo,
            zhi,
        )

    # One oscillator bead per reservoir.
    for i, (cx, cy, cz) in enumerate(centers):
        atoms.append((4 + i, cx, cy, cz))

    # Small Ar reservoirs, excluded from the immediate bead core.
    for i, (cx, cy, cz) in enumerate(centers):
        add_fcc_box(
            atoms,
            1 + i,
            args.ar_lattice,
            i * lx + args.wall_margin,
            (i + 1) * lx - args.wall_margin,
            ylo + args.wall_margin,
            yhi - args.wall_margin,
            zlo + args.wall_margin,
            zhi - args.wall_margin,
            limit=args.ar_per_reservoir,
            exclude_spheres=[(cx, cy, cz, 5.0)],
        )

    box = (xlo, xhi, ylo, yhi, zlo, zhi)
    write_data(Path(args.output), atoms, box, args.bead_mass)

    counts = {i: 0 for i in range(1, 8)}
    for atom_type, *_ in atoms:
        counts[atom_type] += 1
    print(f"Wrote {args.output}")
    print(f"Total atoms: {len(atoms)}")
    for atom_type in range(1, 8):
        print(f"type {atom_type}: {counts[atom_type]}")


if __name__ == "__main__":
    main()
