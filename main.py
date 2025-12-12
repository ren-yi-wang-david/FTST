from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
PY_SRC = ROOT / "src" / "python"
if str(PY_SRC) not in sys.path:
    sys.path.insert(0, str(PY_SRC))

from FTST import ThermalSolver


def main() -> None:
    """Minimal end-to-end FTST workflow for users running python main.py."""
    solver = ThermalSolver(Lx=0.01, Ly=0.01, Lz=0.01, dx=1e-4)
    print(f"[INFO] Grid size = {solver.nx} × {solver.ny} × {solver.nz}")

    solver.set_boundary([25.0, 25.0, 25.0, 25.0, 25.0, 60.0])
    print("[INFO] Boundary set (top face 60 °C, others 25 °C)")

    print("[INFO] Solving ...")
    solver.solve(max_iter=5000, tolerance=1e-8)
    print("[INFO] Solve complete")

    phi = solver.solution_array()
    cx, cy, cz = solver.nx // 2, solver.ny // 2, solver.nz // 2
    print(f"[INFO] Center temperature = {phi[cx, cy, cz]:.4f} °C")

    out_dir = ROOT / "output"
    out_dir.mkdir(exist_ok=True)

    csv_path = out_dir / "result_from_main.csv"
    solver.temp_csv(csv_path.as_posix())
    print(f"[INFO] Temperature field saved → {csv_path}")

    solver.cutplane_plot("z", cz, (out_dir / "heatmap_z_main.png").as_posix())
    print("[INFO] Saved heatmap → heatmap_z_main.png")

    ice_path = ROOT / "tests" / "resources" / "ice.txt"
    if ice_path.exists():
        print(f"[INFO] Running compare() against {ice_path} ...")
        solver.compare(ice_path.as_posix(), csv_path.as_posix())
    else:
        print("[INFO] Icepak reference not found; skipping comparison.")


if __name__ == "__main__":
    main()
