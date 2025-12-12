import numpy as np

from ftst_dense import Dense3DSolver


def make_solver(n, dx=1e-3):
    length = (n - 1) * dx
    solver = Dense3DSolver(length, length, length, dx)
    rhs_shape = (solver.nx, solver.ny, solver.nz)
    solver.set_rhs(np.zeros(rhs_shape, dtype=float))
    return solver


def test_constant_boundary_returns_constant_solution():
    solver = make_solver(6)
    solver.set_boundary([50.0] * 6)
    solver.solve(2000, 1e-6)

    sol = solver.solution_array()
    assert np.allclose(sol, 50.0, atol=1e-3)


def test_gradient_boundary_center_is_between_extremes():
    solver = make_solver(8)
    solver.set_boundary([0.0, 25.0, 50.0, 75.0, 100.0, 125.0])
    solver.solve(4000, 1e-5)

    sol = solver.solution_array()
    center = sol[4, 4, 4]
    assert 0.0 < center < 125.0
    assert np.isfinite(sol).all()
