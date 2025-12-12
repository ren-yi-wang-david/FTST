#pragma once
#include <array>
#include <vector>

#include <ftst/equation.hpp>

using Matrix3DDouble = Matrix3D<double, TrackingAllocator<double>>;

class Dense3DSolver
{
public:
    Dense3DSolver(size_t nx, size_t ny, size_t nz, double dx);
    Dense3DSolver(double lx, double ly, double lz, double dx);

    void setRHS(const Matrix3DDouble &rhs);
    void setBoundary(const std::array<double, 6> &bc);

    void solve(size_t maxIter, double tol);

    const Matrix3DDouble &solution() const { return phi_; }

    size_t nx() const { return nx_; }
    size_t ny() const { return ny_; }
    size_t nz() const { return nz_; }
    double dx() const { return dx_; }

private:
    // --------------------
    // Grid size
    // --------------------
    size_t nx_, ny_, nz_;
    double dx_;
    double omega;
    bool isBoundary(size_t idx) const;
    // --------------------
    // Original 3D matrices (kept for external interface)
    // --------------------
    Matrix3DDouble phi_; // answer (for user)
    Matrix3DDouble rhs_;

    // ==================================================
    // 1D flattened buffers for HPC SOR (internal)
    // ==================================================
    std::vector<double> phi1D_; // flattened phi
    std::vector<double> rhs1D_; // flattened rhs

    // Red/Black lists: store *1D index*, not (i,j,k)
    std::vector<size_t> red_list_;
    std::vector<size_t> black_list_;

    // ==================================================
    // Internal helpers
    // ==================================================
    void buildColoring(); // Build tiled red-black lists (1D index)
    void flattenMemory(); // Allocate & init phi1D_ and rhs1D_
};
