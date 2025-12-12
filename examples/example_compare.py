import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PY_SRC = ROOT / "src" / "python"
if str(PY_SRC) not in sys.path:
    sys.path.insert(0, str(PY_SRC))

from FTST.compare import compare_ice_to_ftst


def main():
    csv_path = ROOT / "output" / "result.csv"
    ice_path = ROOT / "tests" / "resources" / "ice.txt"

    if not csv_path.exists():
        raise FileNotFoundError("Run example_basic.py to generate output/result.csv first.")
    if not ice_path.exists():
        raise FileNotFoundError("Missing tests/resources/ice.txt for comparison.")

    compare_ice_to_ftst(ice_path.as_posix(), csv_path.as_posix(), ROOT / "output")


if __name__ == "__main__":
    main()
