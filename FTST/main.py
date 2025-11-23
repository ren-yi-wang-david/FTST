import os
import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
import FTST

INPUT_DIR = os.path.join(os.path.dirname(__file__), "input")

def run_pipeline():
    """主流程：依序載入幾何、材料、邊界、功率，呼叫求解"""
    FTST.load_geometry(os.path.join(INPUT_DIR, "geometry.json"))
    FTST.load_material(os.path.join(INPUT_DIR, "materials.json"))
    FTST.load_boundary(os.path.join(INPUT_DIR, "boundary.json"))
    FTST.load_power(os.path.join(INPUT_DIR, "power_map.csv"))

    T = np.array(FTST.solve())
    Nx, Ny, Nz =  FTST.load_geometry.Nx, FTST.load_geometry.Ny, FTST.load_geometry.Nz
    T = T.reshape(Nx, Ny, Nz)
    visualize_cutplane(T)

def visualize_cutplane(T):
    Nx, Ny, Nz = T.shape
    print("\n📸 顯示模式選擇：")
    print("  [0] 中心切面 (預設)")
    print("  [1] 手動輸入 Z-index")
    print("  [2] 三向切面")
    mode = input("選擇模式 [0/1/2]: ") or "0"

    if mode == "0":
        z_idx = Nz // 2
        plt.imshow(T[:, :, z_idx], cmap="hot", origin="lower")
        plt.title(f"FTST Cut Z={z_idx}")
        plt.colorbar()
        plt.show()

    elif mode == "1":
        z_idx = int(input(f"輸入 Z-index (0~{Nz-1}): ") or Nz//2)
        plt.imshow(T[:, :, z_idx], cmap="hot", origin="lower")
        plt.title(f"FTST Cut Z={z_idx}")
        plt.colorbar()
        plt.show()

    elif mode == "2":
        fig, axes = plt.subplots(1, 3, figsize=(12, 4))
        axes[0].imshow(T[:, :, Nz//2], cmap="hot", origin="lower")
        axes[0].set_title(f"XY (Z={Nz//2})")
        axes[1].imshow(T[:, Ny//2, :], cmap="hot", origin="lower")
        axes[1].set_title(f"XZ (Y={Ny//2})")
        axes[2].imshow(T[Nx//2, :, :], cmap="hot", origin="lower")
        axes[2].set_title(f"YZ (X={Nx//2})")
        for ax in axes: ax.axis("off")
        plt.suptitle("FTST Multi-View")
        plt.tight_layout()
        plt.show()

if __name__ == "__main__":
    run_pipeline()
