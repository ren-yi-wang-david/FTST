import numpy as np

from ftst_dense import Dense3DSolver


def test_boundary_faces_match_requested_temperatures():
    solver = Dense3DSolver(0.003, 0.003, 0.003, 1e-3)
    solver.set_rhs(np.zeros((solver.nx, solver.ny, solver.nz), dtype=float))

    boundaries = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
    solver.set_boundary(boundaries)
    solver.solve(1500, 1e-6)

    sol = solver.solution_array()

    def assert_face(face, expected):
        view = face
        if view.shape[0] > 2:
            view = view[1:-1, :]
        if view.shape[1] > 2:
            view = view[:, 1:-1]
        assert np.allclose(view, expected, atol=1e-2)

    assert_face(sol[0, :, :], boundaries[0])
    assert_face(sol[-1, :, :], boundaries[1])
    assert_face(sol[:, 0, :], boundaries[2])
    assert_face(sol[:, -1, :], boundaries[3])
    assert_face(sol[:, :, 0], boundaries[4])
    assert_face(sol[:, :, -1], boundaries[5])
