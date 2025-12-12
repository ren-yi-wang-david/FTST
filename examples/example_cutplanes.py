from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
CSV_PATH = ROOT / "output" / "result.csv"

if not CSV_PATH.exists():
    raise FileNotFoundError(
        f"{CSV_PATH} not found. Generate it first via examples/example_basic.py."
    )

# ==========================================================
# 讀取 temperature_xyz.csv（格式：x_m, y_m, z_m, phi）
# ==========================================================
raw = np.loadtxt(CSV_PATH, delimiter=",", skiprows=1)

xs_all  = raw[:, 0]
ys_all  = raw[:, 1]
zs_all  = raw[:, 2]
temps_all = raw[:, 3]


# ----------------------------------------------------------
# 產生 unique 座標索引表（自動）
# ----------------------------------------------------------
ux = np.unique(xs_all)
uy = np.unique(ys_all)
uz = np.unique(zs_all)

# 座標 → index 映射
x_id = {v: i for i, v in enumerate(ux)}
y_id = {v: i for i, v in enumerate(uy)}
z_id = {v: i for i, v in enumerate(uz)}

NX = len(ux)
NY = len(uy)
NZ = len(uz)

print(f"[INFO] NX={NX}, NY={NY}, NZ={NZ}")
print(f"[INFO] dx = {ux[1] - ux[0]} m")

# ==========================================================
# 設切面
# mode: 'x', 'y', 'z'
# target_m: 切面位置（公尺）
# ==========================================================
mode = 'z'
target_m = 0.005  # 5 mm

# 找到最接近 target_m 的座標值
def nearest(arr, v):
    idx = np.argmin(np.abs(arr - v))
    return idx, arr[idx]

if mode == 'z':
    level_index, target = nearest(uz, target_m)
elif mode == 'y':
    level_index, target = nearest(uy, target_m)
elif mode == 'x':
    level_index, target = nearest(ux, target_m)

print(f"[INFO] mode={mode}, target={target_m} m, actual={target} m")

# ==========================================================
# 建立切面矩陣
# ==========================================================
if mode == 'z':
    plane = np.zeros((NY, NX))
    mask = (zs_all == target)

    for x, y, t in zip(xs_all[mask], ys_all[mask], temps_all[mask]):
        plane[y_id[y], x_id[x]] = t

    xlabel, ylabel = "x (m)", "y (m)"
    extent = [ux.min(), ux.max(), uy.min(), uy.max()]
    title = f"Temperature Cut-plane at z={target:.4f} m"

elif mode == 'y':
    plane = np.zeros((NZ, NX))
    mask = (ys_all == target)

    for x, z, t in zip(xs_all[mask], zs_all[mask], temps_all[mask]):
        plane[z_id[z], x_id[x]] = t

    xlabel, ylabel = "x (m)", "z (m)"
    extent = [ux.min(), ux.max(), uz.min(), uz.max()]
    title = f"Temperature Cut-plane at y={target:.4f} m"

elif mode == 'x':
    plane = np.zeros((NY, NZ))
    mask = (xs_all == target)

    for y, z, t in zip(ys_all[mask], zs_all[mask], temps_all[mask]):
        plane[y_id[y], z_id[z]] = t

    xlabel, ylabel = "z (m)", "y (m)"
    extent = [uz.min(), uz.max(), uy.min(), uy.max()]
    title = f"Temperature Cut-plane at x={target:.4f} m"

# ==========================================================
# Plot
# ==========================================================
plt.figure(figsize=(7, 6))
plt.imshow(
    plane,
    origin='lower',
    extent=extent,
    aspect='auto',
    cmap="inferno"
)
plt.colorbar(label="Temperature (°C)")
plt.title(title)
plt.xlabel(xlabel)
plt.ylabel(ylabel)

plt.tight_layout()
plt.show()
