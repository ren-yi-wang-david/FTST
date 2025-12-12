import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PY_SRC = ROOT / "src" / "python"
if str(PY_SRC) not in sys.path:
    sys.path.insert(0, str(PY_SRC))

from FTST import ThermalSolver


def main():

    # ============================================================
    # 1. 建立 Solver
    # ============================================================
    solver = ThermalSolver(Lx=0.01, Ly=0.01, Lz=0.01, dx=1e-4)

    print(f"Grid size: {solver.nx} × {solver.ny} × {solver.nz}")
    print(f"dx = {solver.dx}")

    # ============================================================
    # 2. 設定六面邊界條件
    #    (x-, x+, y-, y+, z-, z+)
    # ============================================================
    solver.set_boundary(25.0, 25.0, 25.0, 25.0, 25.0, 60.0)
    print("Boundary set: 25/25/25/25/25/60 °C")

    cx, cy, cz = solver.nx // 2, solver.ny // 2, solver.nz // 2

    # ============================================================
    # 4. 求解 (Gauss–Seidel + Red–Black + AVX2 + OpenMP)
    # ============================================================
    print("Solving...")
    solver.solve(max_iter=5000, tolerance=1e-8)
    print("Solve finished!")

    # ============================================================
    # 5. 匯出 CSV
    # ============================================================
    out_dir = ROOT / "output"
    out_dir.mkdir(exist_ok=True)

    csv_path = out_dir / "result.csv"
    solver.temp_csv(csv_path.as_posix())
    print(f"Temperature field saved → {csv_path}")

    # ============================================================
    # 6. 匯出切面圖：z-plane
    # ============================================================
    solver.cutplane_plot("z", cz, (out_dir / "heatmap_z.png").as_posix())
    solver.cutplane_plot("x", cx, (out_dir / "heatmap_x.png").as_posix())
    solver.cutplane_plot("y", cy, (out_dir / "heatmap_y.png").as_posix())
    print("Cut-plane heatmaps saved in output/")

    # ============================================================
    # 7. 取得 numpy array (可做後處理)
    # ============================================================
    phi = solver.solution_array()
    print("Center temperature =", phi[cx, cy, cz])

    # ============================================================
    # 8. 與 Icepak 結果比較 (若提供 ICE.csv)
    # ============================================================
    ice_path = ROOT / "tests" / "resources" / "ice.txt"
    if ice_path.exists():
        print("Running compare() against ice.txt ...")
        result = solver.compare(ice_path.as_posix(), csv_path.as_posix())
        print("Compare result:", result)
    else:
        print("No ice.txt found → skipping compare().")


if __name__ == "__main__":
    main()
