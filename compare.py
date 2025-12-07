import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata


# ==========================================================
# Icepak：2D cut-plane（單位 = m）
# Icepak 檔案格式通常是 whitespace 分隔
# ==========================================================
def load_icepak(path):
    df = pd.read_csv(
        path,
        delim_whitespace=True,
        skiprows=3,
        names=["x_m", "y_m", "z_m", "T_ice"]
    )
    return df


# ==========================================================
# FTST：3D temperature_xyz.csv（單位 = m）
# 你的 C++ 格式為：x_m, y_m, z_m, phi
# ==========================================================
def load_ftst(path):
    df = pd.read_csv(path)

    # 統一欄位名稱
    df.rename(columns={
        "x_m": "x",
        "y_m": "y",
        "z_m": "z",
        "phi": "T_ftst"
    }, inplace=True)

    return df


# ==========================================================
# 找出 FTST 中最接近 Icepak z-plane
# ==========================================================
def get_matching_ftst_slice(ftst, ice):
    target_z = ice["z_m"].iloc[0]  # Icepak 固定 z-plane
    uz = np.unique(ftst["z"])

    idx = np.argmin(np.abs(uz - target_z))
    matched_z = uz[idx]

    print(f"🔎 Ice z = {target_z} m → 使用 FTST z = {matched_z} m")

    ftst_slice = ftst[np.abs(ftst["z"] - matched_z) < 1e-12].copy()
    return ftst_slice, matched_z


# ==========================================================
# 插值 FTST → Icepak 座標
# ==========================================================
def interpolate_ftst_to_ice(ice, ftst_slice):
    pts = ftst_slice[["x", "y"]].values
    vals = ftst_slice["T_ftst"].values

    target_pts = ice[["x_m", "y_m"]].values

    interp = griddata(pts, vals, target_pts, method="linear")

    return interp


# ==========================================================
# 畫 heatmap
# ==========================================================
def plot_map(ice, field, title):
    pivot = ice.pivot(index="y_m", columns="x_m", values=field)

    plt.figure(figsize=(7, 6))
    plt.imshow(
        pivot.values,
        origin="lower",
        extent=[
            pivot.columns.min(), pivot.columns.max(),
            pivot.index.min(), pivot.index.max()
        ],
        cmap="inferno",
        aspect="auto"
    )
    plt.colorbar(label=field)
    plt.title(title)
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.tight_layout()
    plt.show()


# ==========================================================
# 主比較流程（Icepak 視為標準）
# ==========================================================
def compare_ice_to_ftst(ice_file, ftst_file):
    ice = load_icepak(ice_file)
    ftst = load_ftst(ftst_file)

    print(f"📘 Icepak points = {len(ice)}")
    print(f"📗 FTST points = {len(ftst)}")

    ftst_slice, matched_z = get_matching_ftst_slice(ftst, ice)

    interp_ftst = interpolate_ftst_to_ice(ice, ftst_slice)
    ice["T_ftst_interp"] = interp_ftst

    ice["abs_err"] = np.abs(ice["T_ice"] - ice["T_ftst_interp"])
    ice["pct_err"] = ice["abs_err"] / (np.abs(ice["T_ice"]) + 1e-12) * 100

    print("\n============================")
    print("📊 Error Summary (Ice as standard)")
    print("============================")
    print(f"Max abs error  : {ice['abs_err'].max():.6f} °C")
    print(f"Mean abs error : {ice['abs_err'].mean():.6f} °C")
    print(f"Max % error    : {ice['pct_err'].max():.3f} %")
    print(f"Mean % error   : {ice['pct_err'].mean():.3f} %")

    plot_map(ice, "abs_err", "Absolute Error Heatmap (Ice as standard)")
    plot_map(ice, "pct_err", "Percentage Error Heatmap (Ice as standard)")

    ice.to_csv("difference_map_IceAsStandard.csv", index=False)
    print("\n📁 已輸出 difference_map_IceAsStandard.csv")

    return ice


if __name__ == "__main__":
    compare_ice_to_ftst("ice.txt", "phi_output.csv")
