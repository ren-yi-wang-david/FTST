from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

from ftst_dense import Dense3DSolver

from .compare import compare_ice_to_ftst
from .utils.io import write_temperature_csv
from .utils.plotting import save_heatmap


class ThermalSolver:
    """
    Thin convenience wrapper around the pybind11 Dense3DSolver.
    Provides a friendlier API for Python scripts.
    """

    def __init__(self, *, Lx: float, Ly: float, Lz: float, dx: float):
        self._solver = Dense3DSolver(Lx, Ly, Lz, dx)
        self._set_zero_rhs()

    # -----------------------------
    # Basic properties
    # -----------------------------
    @property
    def nx(self) -> int:
        return self._solver.nx

    @property
    def ny(self) -> int:
        return self._solver.ny

    @property
    def nz(self) -> int:
        return self._solver.nz

    @property
    def dx(self) -> float:
        return self._solver.dx

    # -----------------------------
    # Solver configuration
    # -----------------------------
    def _set_zero_rhs(self) -> None:
        shape = (self.nx, self.ny, self.nz)
        self._solver.set_rhs(np.zeros(shape, dtype=float))

    def set_boundary(self, *values: float) -> None:
        """
        Accept six boundary values (x-, x+, y-, y+, z-, z+).
        Also accepts a single iterable of length six.
        """
        if len(values) == 1 and isinstance(values[0], Iterable):
            values = tuple(values[0])  # type: ignore[assignment]
        if len(values) != 6:
            raise ValueError("Expected six boundary values: x-, x+, y-, y+, z-, z+")
        bc = [float(v) for v in values]
        self._solver.set_boundary(bc)

    def solve(self, max_iter: int = 10000, tolerance: float = 1e-9) -> None:
        self._solver.solve(int(max_iter), float(tolerance))

    def solution_array(self) -> np.ndarray:
        return self._solver.solution_array()

    # -----------------------------
    # Outputs
    # -----------------------------
    def temp_csv(self, filename: str) -> None:
        """
        Dump the temperature field to CSV with columns x_m,y_m,z_m,phi.
        """
        write_temperature_csv(
            self.solution_array(),
            self.dx,
            Path(filename),
        )

    def cutplane_plot(self, axis: str, index: int, filename: str | None = None, cmap: str = "inferno") -> None:
        """
        Save a heatmap of a temperature cut-plane.
        axis: 'x', 'y', or 'z'
        index: slice index along the chosen axis
        filename: optional; if None, auto-generate heatmap_<axis>_<index>.png
        """
        axis_clean = axis.lower()
        if axis_clean not in {"x", "y", "z"}:
            raise ValueError("axis must be 'x', 'y', or 'z'")
        idx = int(index)

        # Auto filename
        if filename is None:
            filename = f"heatmap_{axis_clean}_{idx}.png"

        path = Path(filename)
        if path.parent and not path.parent.exists():
            path.parent.mkdir(parents=True, exist_ok=True)

        plane, extent, xlabel, ylabel, title = self._make_slice(axis_clean, idx)

        save_heatmap(plane, extent, xlabel, ylabel, title, path, cmap)

        axis_label = {"x": "x-cutplane", "y": "y-cutplane", "z": "z-cutplane"}[axis_clean]
        print(f"[INFO] Saved {axis_label} (index={idx}) → {path.as_posix()}")

    def compare(self, ice_file: str, ftst_csv: str | None = None):
        """Run the compare.py pipeline against the current solution."""
        ice_path = Path(ice_file)
        if not ice_path.exists():
            raise FileNotFoundError(f"Icepak reference file not found: {ice_path}")

        csv_path = Path(ftst_csv) if ftst_csv else Path("output/result.csv")
        # Always write the latest field to ensure we compare against the current solve.
        self.temp_csv(str(csv_path))

        return compare_ice_to_ftst(str(ice_path), str(csv_path))

    # -----------------------------
    # Helpers
    # -----------------------------
    def _axis_span(self, count: int) -> float:
        return 0.0 if count <= 1 else (count - 1) * self.dx

    def _make_slice(self, axis: str, index: int):
        phi = self.solution_array()
        if axis == "z":
            if not 0 <= index < self.nz:
                raise IndexError("z index out of range")
            plane = phi[:, :, index].T
            extent = (0.0, self._axis_span(self.nx), 0.0, self._axis_span(self.ny))
            xlabel, ylabel = "x (m)", "y (m)"
            title = f"Temperature Cut-plane at z={index * self.dx} m"
        elif axis == "y":
            if not 0 <= index < self.ny:
                raise IndexError("y index out of range")
            plane = phi[:, index, :].T
            extent = (0.0, self._axis_span(self.nx), 0.0, self._axis_span(self.nz))
            xlabel, ylabel = "x (m)", "z (m)"
            title = f"Temperature Cut-plane at y={index * self.dx} m"
        else:  # axis == "x"
            if not 0 <= index < self.nx:
                raise IndexError("x index out of range")
            plane = phi[index, :, :]
            extent = (0.0, self._axis_span(self.nz), 0.0, self._axis_span(self.ny))
            xlabel, ylabel = "z (m)", "y (m)"
            title = f"Temperature Cut-plane at x={index * self.dx} m"

        return plane, extent, xlabel, ylabel, title


__all__: Sequence[str] = ("ThermalSolver",)
