import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PY_SRC = ROOT / "src" / "python"
if str(PY_SRC) not in sys.path:
    sys.path.insert(0, str(PY_SRC))

from FTST.compare import compare_ice_to_ftst


def _write_ice_file(path: Path, entries: list[tuple[float, float, float, float]]) -> None:
    with path.open("w", encoding="utf-8") as fout:
        fout.write("1\n")
        fout.write("units m C\n")
        fout.write(f"{len(entries)}\n")
        for x, y, z, t in entries:
            fout.write(f"{x} {y} {z} {t}\n")


def _write_ftst_csv(path: Path, entries: list[tuple[float, float, float, float]]) -> None:
    with path.open("w", encoding="utf-8") as fout:
        fout.write("x_m,y_m,z_m,phi\n")
        for x, y, z, t in entries:
            fout.write(f"{x},{y},{z},{t}\n")


def test_compare_detects_zero_error(tmp_path):
    common_points = [
        (0.0, 0.0, 0.01, 25.0),
        (0.0, 0.001, 0.01, 30.0),
        (0.001, 0.0, 0.01, 27.5),
        (0.001, 0.001, 0.01, 29.5),
    ]
    ice_file = tmp_path / "ice.txt"
    ftst_file = tmp_path / "field.csv"

    _write_ice_file(ice_file, common_points)
    _write_ftst_csv(ftst_file, common_points)

    df = compare_ice_to_ftst(ice_file.as_posix(), ftst_file.as_posix(), tmp_path)
    assert np.isclose(df["abs_err"].max(), 0.0)
    assert np.isclose(df["pct_err"].max(), 0.0)
