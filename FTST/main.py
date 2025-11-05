import numpy as np
import FTST

Nx, Ny = 10, 10
dx = dy = 1.0
k = 1.0

power = np.zeros(Nx * Ny)
bc = np.zeros(Nx * Ny)
power[(Nx//2)*Ny + (Ny//2)] = 10.0  # 中央熱源

T = FTST.solve(Nx, Ny, dx, dy, k, power, bc)
print("Temperature field:")
print(T.reshape(Nx, Ny))
