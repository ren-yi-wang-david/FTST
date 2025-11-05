#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>

namespace ftst {

class ThermalAnalysis {
public:
    ThermalAnalysis(int nx, int ny, double dx, double dy, double k)
        : Nx(nx), Ny(ny), dx(dx), dy(dy), k(k), 
          power(nx * ny, 0.0), bc(nx * ny, 0.0), T(nx * ny, 0.0) {}

    void set_power(const std::vector<double>& p) {
        if (p.size() != power.size()) throw std::runtime_error("Power size mismatch");
        power = p;
    }

    void set_boundary(const std::vector<double>& b) {
        if (b.size() != bc.size()) throw std::runtime_error("Boundary size mismatch");
        bc = b;
    }

    void solve_steady(int iterations = 5000) {
        // 五點格式 Jacobi 迭代求解
        for (int it = 0; it < iterations; ++it) {
            for (int i = 1; i < Nx - 1; ++i) {
                for (int j = 1; j < Ny - 1; ++j) {
                    int idx = i * Ny + j;
                    T[idx] = 0.25 * (T[(i+1)*Ny+j] + T[(i-1)*Ny+j] +
                                     T[i*Ny+(j+1)] + T[i*Ny+(j-1)]) +
                              power[idx] * dx * dy / (4 * k);
                }
            }
        }
    }

    const std::vector<double>& get_temperature() const { return T; }

private:
    int Nx, Ny;
    double dx, dy, k;
    std::vector<double> power, bc, T;
};

} // namespace ftst
