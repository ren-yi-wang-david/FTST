from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np


def _ensure_parent(path: Path) -> Path:
    if path.parent and not path.parent.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
    return path


def write_temperature_csv(
    phi: np.ndarray,
    dx: float,
    destination: Path | str,
) -> None:
    """Persist the temperature volume as a CSV with x,y,z indices scaled by dx."""
    nx, ny, nz = phi.shape
    path = _ensure_parent(Path(destination))

    with path.open("w", encoding="utf-8") as fout:
        fout.write("x_m,y_m,z_m,phi\n")
        for i in range(nx):
            x = i * dx
            for j in range(ny):
                y = j * dx
                for k in range(nz):
                    z = k * dx
                    fout.write(f"{x},{y},{z},{phi[i, j, k]}\n")


__all__: Iterable[str] = ("write_temperature_csv",)
