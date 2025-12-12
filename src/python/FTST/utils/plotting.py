from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np


def save_heatmap(
    plane: np.ndarray,
    extent: Sequence[float],
    xlabel: str,
    ylabel: str,
    title: str,
    destination: Path | str,
    cmap: str = "inferno",
) -> None:
    """Persist a heatmap image for the provided plane data."""
    path = Path(destination)
    if path.parent and not path.parent.exists():
        path.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(7, 6))
    plt.imshow(
        plane,
        origin="lower",
        extent=extent,
        aspect="auto",
        cmap=cmap,
    )
    plt.colorbar(label="Temperature (°C)")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()


__all__: Iterable[str] = ("save_heatmap",)
