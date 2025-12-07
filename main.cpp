#include <iostream>
#include <array>
#include <chrono>
#include <fstream>
#include "solver.hpp"
#include "memory_tracker.hpp"

// ======================================================================
//  儲存 CSV
// ======================================================================
void saveCSV(const Matrix3D<double> &phi,
             size_t Nx, size_t Ny, size_t Nz,
             double dx_m,
             const std::string &filename)
{
    auto t0 = std::chrono::high_resolution_clock::now();

    std::ofstream fout(filename);
    if (!fout)
    {
        std::cerr << "Cannot open file: " << filename << "\n";
        return;
    }

    fout << "x_m,y_m,z_m,phi\n";

    for (size_t i = 0; i < Nx; i++)
        for (size_t j = 0; j < Ny; j++)
            for (size_t k = 0; k < Nz; k++)
            {
                double x = i * dx_m;
                double y = j * dx_m;
                double z = k * dx_m;

                fout << x << "," << y << "," << z << "," << phi(i, j, k) << "\n";
            }

    fout.close();

    auto t1 = std::chrono::high_resolution_clock::now();
    long long ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    std::cout << "CSV saved in " << ms << " ms\n";
}

// ======================================================================
//  主程式
// ======================================================================
int main()
{
    double Lx = 0.01; // 10 mm in meters
    double Ly = 0.01; // 10 mm in meters
    double Lz = 0.01; // 10 mm in meters

    double dx = 1e-4; // 0.1 mm in meters

    std::cout << "Before allocation:\n";
    MemoryTracker::report();

    {
        auto t_alloc0 = std::chrono::high_resolution_clock::now();

        Dense3DSolver solver(Lx, Ly, Lz, dx);
        size_t Nx = solver.nx();
        size_t Ny = solver.ny();
        size_t Nz = solver.nz();

        auto t_alloc1 = std::chrono::high_resolution_clock::now();
        long long alloc_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                                 t_alloc1 - t_alloc0)
                                 .count();
        std::cout << "Solver allocation time = " << alloc_ms << " ms\n";

        MemoryTracker::report();

        // ==============================================================
        // 建立 RHS
        // ==============================================================
        auto t_rhs0 = std::chrono::high_resolution_clock::now();

        Matrix3D<double> rhs(Nx, Ny, Nz);
        for (size_t i = 0; i < Nx; i++)
            for (size_t j = 0; j < Ny; j++)
                for (size_t k = 0; k < Nz; k++)
                    rhs(i, j, k) = 0.0;

        solver.setRHS(rhs);

        auto t_rhs1 = std::chrono::high_resolution_clock::now();
        long long rhs_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                               t_rhs1 - t_rhs0)
                               .count();
        std::cout << "RHS setup time = " << rhs_ms << " ms\n";

        MemoryTracker::report();

        // ==============================================================
        // 邊界條件
        // ==============================================================
        std::array<double, 6> bc = {25.0, 25.0, 25.0, 25.0, 25.0, 60.0};

        auto t_bc0 = std::chrono::high_resolution_clock::now();
        solver.setBoundary(bc);
        auto t_bc1 = std::chrono::high_resolution_clock::now();

        long long bc_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                              t_bc1 - t_bc0)
                              .count();
        std::cout << "Boundary set time = " << bc_ms << " ms\n";

        MemoryTracker::report();

        // ==============================================================
        // 求解 (Gauss-Seidel / LU fallback 依你的 solver 實作而定)
        // ==============================================================

        auto t_solve0 = std::chrono::high_resolution_clock::now();
        solver.solve(100000, 1e-10);
        auto t_solve1 = std::chrono::high_resolution_clock::now();

        long long solve_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                                 t_solve1 - t_solve0)
                                 .count();
        std::cout << "Solve time = " << solve_ms << " ms\n";

        MemoryTracker::report();

        // ==============================================================
        // 中心點結果
        // ==============================================================
        const auto &phi = solver.solution();
        std::cout << "phi(center) = "
                  << phi(Nx / 2, Ny / 2, Nz / 2)
                  << "\n";

        // ==============================================================
        // 輸出 CSV
        // ==============================================================
        saveCSV(phi, Nx, Ny, Nz, dx, "phi_output.csv");

        std::cout << "\nInside scope memory:\n";
        MemoryTracker::report();
    }

    std::cout << "\nAfter scope (should be freed):\n";
    MemoryTracker::report();

    return 0;
}
