#pragma once
#include "equation.hpp"

using Matrix3DDouble = Matrix3D<double, TrackingAllocator<double>>;

class Dense3DSolver
{
public:
    Dense3DSolver(size_t nx, size_t ny, size_t nz, double dx);
    Dense3DSolver(double lx, double ly, double lz, double dx);

    // 設定 RHS（Poisson 方程中的來源項）
    void setRHS(const Matrix3DDouble &rhs);

    // 設定 Dirichlet 邊界條件
    void setBoundary(const std::array<double, 6> &bc);

    // Gauss-Seidel 求解 Laplace/Poisson
    void solve(size_t maxIter, double tol);

    // 取得解
    const Matrix3DDouble &solution() const { return phi_; }

    size_t nx() const { return nx_; }
    size_t ny() const { return ny_; }
    size_t nz() const { return nz_; }
    double dx() const { return dx_; }

private:
    size_t nx_, ny_, nz_;
    double dx_;
    double omega ;
    Matrix3DDouble phi_; // 解
    Matrix3DDouble rhs_; // 來源項（Poisson 方程）

};
