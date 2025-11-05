import numpy as np
import pytest
import FTST


def test_import():
    """確認模組能正常匯入"""
    assert hasattr(FTST, "solve"), "FTST 模組中找不到 solve() 函式"


def test_solver_shape():
    """確認輸出結果大小正確"""
    Nx, Ny = 6, 6
    dx = dy = 1.0
    k = 1.0
    power = np.zeros(Nx * Ny)
    bc = np.zeros(Nx * Ny)
    T = FTST.solve(Nx, Ny, dx, dy, k, power, bc)
    assert isinstance(T, np.ndarray)
    assert T.size == Nx * Ny


def test_center_hotspot():
    """確認中央熱源產生溫度升高"""
    Nx, Ny = 8, 8
    dx = dy = 1.0
    k = 1.0
    power = np.zeros(Nx * Ny)
    bc = np.zeros(Nx * Ny)

    # 在中央放一個熱源
    center = (Nx // 2) * Ny + (Ny // 2)
    power[center] = 10.0

    T = FTST.solve(Nx, Ny, dx, dy, k, power, bc)
    T2D = T.reshape(Nx, Ny)

    # 驗證中心比邊界熱
    center_T = T2D[Nx // 2, Ny // 2]
    edge_T = np.mean(T2D[0, :] + T2D[-1, :] + T2D[:, 0] + T2D[:, -1])
    assert center_T > edge_T, "中心溫度應高於邊界"


def test_symmetric_result():
    """確認溫度分佈對稱"""
    Nx, Ny = 8, 8
    dx = dy = 1.0
    k = 1.0
    power = np.zeros(Nx * Ny)
    bc = np.zeros(Nx * Ny)
    power[(Nx // 2) * Ny + (Ny // 2)] = 10.0

    T = FTST.solve(Nx, Ny, dx, dy, k, power, bc)
    T2D = T.reshape(Nx, Ny)

    # 檢查對稱性
    diff = np.abs(T2D - np.flipud(np.fliplr(T2D)))
    assert np.max(diff) < 1e-6, "溫度分佈應該是對稱的"
