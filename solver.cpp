
#include "solver.hpp"
#include <cmath>
#include <array>
#include <omp.h>

namespace
{
size_t to_steps(double length, double spacing)
{
    return static_cast<size_t>(std::llround(length / spacing)) + 1;
}
}
Dense3DSolver::Dense3DSolver(size_t nx, size_t ny, size_t nz, double dx)
    : nx_(nx), ny_(ny), nz_(nz), dx_(dx),
      phi_(nx, ny, nz),
      rhs_(nx, ny, nz)
{
    // 初始值
    for (size_t i = 0; i < nx_; i++)
        for (size_t j = 0; j < ny_; j++)
            for (size_t k = 0; k < nz_; k++)
                phi_(i, j, k) = 0.0;

    // -------------------------
    // 計算最佳 ω
    // -------------------------
    double rx = std::cos(M_PI / nx_);
    double ry = std::cos(M_PI / ny_);
    double rz = std::cos(M_PI / nz_);
    double rho = (rx + ry + rz) / 3.0;

    omega = 2.0 / (1.0 + std::sqrt(1.0 - rho * rho));
}

Dense3DSolver::Dense3DSolver(double lx, double ly, double lz, double dx)
    : Dense3DSolver(to_steps(lx, dx), to_steps(ly, dx), to_steps(lz, dx), dx)
{
}

void Dense3DSolver::setRHS(const Matrix3DDouble &rhs)
{
    for (size_t i = 0; i < nx_; i++)
        for (size_t j = 0; j < ny_; j++)
            for (size_t k = 0; k < nz_; k++)
                rhs_(i, j, k) = rhs(i, j, k);
}

void Dense3DSolver::setBoundary(const std::array<double, 6> &bc)
{

    for (size_t j = 0; j < ny_; j++)
        for (size_t k = 0; k < nz_; k++)
        {
            phi_(0, j, k) = bc[0];
            phi_(nx_ - 1, j, k) = bc[1];
        }

    for (size_t i = 0; i < nx_; i++)
        for (size_t k = 0; k < nz_; k++)
        {
            phi_(i, 0, k) = bc[2];
            phi_(i, ny_ - 1, k) = bc[3];
        }

    for (size_t i = 0; i < nx_; i++)
        for (size_t j = 0; j < ny_; j++)
        {
            phi_(i, j, 0) = bc[4];
            phi_(i, j, nz_ - 1) = bc[5];
        }
}
void Dense3DSolver::solve(size_t maxIter, double tol)
{
    const double dx2 = dx_ * dx_;
    const double w = omega;

    // Tiling block size
    const int Ti = 16, Tj = 16, Tk = 16;

    for (size_t iter = 0; iter < maxIter; ++iter)
    {
        double maxDiff = 0.0;

// ============================
//           RED SWEEP
// ============================
#pragma omp parallel for collapse(3) reduction(max : maxDiff)
        for (int ii = 1; ii < (int)nx_ - 1; ii += Ti)
            for (int jj = 1; jj < (int)ny_ - 1; jj += Tj)
                for (int kk = 1; kk < (int)nz_ - 1; kk += Tk)
                {
                    for (int i = ii; i < std::min(ii + Ti, (int)nx_ - 1); i++)
                        for (int j = jj; j < std::min(jj + Tj, (int)ny_ - 1); j++)
                            for (int k = kk; k < std::min(kk + Tk, (int)nz_ - 1); k++)
                            {
                                if ((i + j + k) & 1)
                                    continue; // red points

                                const double oldVal = phi_(i, j, k);

                                const double gs =
                                    (phi_(i + 1, j, k) + phi_(i - 1, j, k) +
                                     phi_(i, j + 1, k) + phi_(i, j - 1, k) +
                                     phi_(i, j, k + 1) + phi_(i, j, k - 1) -
                                     dx2 * rhs_(i, j, k)) /
                                    6.0;

                                const double newVal = oldVal + w * (gs - oldVal);
                                phi_(i, j, k) = newVal;

                                maxDiff = std::max(maxDiff, std::abs(newVal - oldVal));
                            }
                }

// ============================
//           BLACK SWEEP
// ============================
#pragma omp parallel for collapse(3) reduction(max : maxDiff)
        for (int ii = 1; ii < (int)nx_ - 1; ii += Ti)
            for (int jj = 1; jj < (int)ny_ - 1; jj += Tj)
                for (int kk = 1; kk < (int)nz_ - 1; kk += Tk)
                {
                    for (int i = ii; i < std::min(ii + Ti, (int)nx_ - 1); i++)
                        for (int j = jj; j < std::min(jj + Tj, (int)ny_ - 1); j++)
                            for (int k = kk; k < std::min(kk + Tk, (int)nz_ - 1); k++)
                            {
                                if (!((i + j + k) & 1))
                                    continue; // black points

                                const double oldVal = phi_(i, j, k);

                                const double gs =
                                    (phi_(i + 1, j, k) + phi_(i - 1, j, k) +
                                     phi_(i, j + 1, k) + phi_(i, j - 1, k) +
                                     phi_(i, j, k + 1) + phi_(i, j, k - 1) -
                                     dx2 * rhs_(i, j, k)) /
                                    6.0;

                                const double newVal = oldVal + w * (gs - oldVal);
                                phi_(i, j, k) = newVal;

                                maxDiff = std::max(maxDiff, std::abs(newVal - oldVal));
                            }
                }

        if (maxDiff < tol)
            return;
    }
}
