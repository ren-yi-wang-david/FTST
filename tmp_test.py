import os
import sys

import numpy as np

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
if CURRENT_DIR not in sys.path:
    sys.path.insert(0, CURRENT_DIR)

import ftst_dense


def main():
    Lx = 0.01  # 20 mm
    Ly = 0.01  # 20 mm
    Lz = 0.01  # 10 mm
    dx = 1e-4  # 0.2 mm grid spacing

    solver = ftst_dense.Dense3DSolver(Lx, Ly, Lz, dx)
    nx, ny, nz = solver.nx, solver.ny, solver.nz

    rhs = np.zeros((nx, ny, nz), dtype=float)
    solver.set_rhs(rhs)

    boundary = [25.0, 25.0, 25.0, 20.0, 25.0, 60.0]
    solver.set_boundary(boundary)

    solver.solve(10000, 1e-9)

    phi = solver.solution_array()
    print("solution shape:", phi.shape)

    ftst_dense.plot_temperature_slice(solver, axis="z", position=Lz / 2)


if __name__ == "__main__":
    main()
