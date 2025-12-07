#include <gtest/gtest.h>

#include "../solver.hpp"
#include "../memory_tracker.hpp"

#include <cmath>
#include <array>
#include <unordered_set>
#include <vector>
#include <chrono>

//
// ============================================================
// Helper functions
// ============================================================
namespace
{
    Matrix3D<double> make_zero_rhs(size_t nx, size_t ny, size_t nz)
    {
        Matrix3D<double> rhs(nx, ny, nz);
        for (size_t i = 0; i < nx; ++i)
            for (size_t j = 0; j < ny; ++j)
                for (size_t k = 0; k < nz; ++k)
                    rhs(i, j, k) = 0.0;
        return rhs;
    }
} // namespace

//
// ============================================================
// Matrix3D tests
// ============================================================
TEST(Matrix3DTest, StoresAndRetrievesValues)
{
    Matrix3D<int> matrix(2, 3, 4);
    matrix(1, 2, 3) = 42;
    EXPECT_EQ(matrix(1, 2, 3), 42);
}

TEST(Matrix3DTest, CopyConstructorClonesValues)
{
    Matrix3D<double> matrix(2, 2, 2);
    matrix(1, 1, 1) = 12.5;

    Matrix3D<double> copy(matrix);
    EXPECT_DOUBLE_EQ(copy(1, 1, 1), 12.5);
}

//
// ============================================================
// Dense3DSolver baseline tests
// ============================================================
TEST(Dense3DSolverTest, ConstantBoundaryProducesConstantSolution)
{
    const size_t n = 6;
    const double dx = 1.0;
    Dense3DSolver solver((n - 1) * dx, (n - 1) * dx, (n - 1) * dx, dx);
    ASSERT_EQ(solver.nx(), n);

    auto rhs = make_zero_rhs(solver.nx(), solver.ny(), solver.nz());
    solver.setRHS(rhs);

    std::array<double, 6> bc = {75.0, 75.0, 75.0, 75.0, 75.0, 75.0};
    solver.setBoundary(bc);

    solver.solve(2000, 1e-6);

    const auto &phi = solver.solution();
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            for (size_t k = 0; k < n; ++k)
                EXPECT_NEAR(phi(i, j, k), 75.0, 1e-3);
}

TEST(Dense3DSolverTest, ZeroBoundaryStaysNearZero)
{
    const size_t n = 5;
    const double dx = 1.0;
    Dense3DSolver solver((n - 1) * dx, (n - 1) * dx, (n - 1) * dx, dx);
    ASSERT_EQ(solver.nx(), n);

    auto rhs = make_zero_rhs(solver.nx(), solver.ny(), solver.nz());
    solver.setRHS(rhs);
    std::array<double, 6> bc = {0, 0, 0, 0, 0, 0};
    solver.setBoundary(bc);

    solver.solve(500, 1e-6);

    const auto &phi = solver.solution();
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            for (size_t k = 0; k < n; ++k)
                EXPECT_NEAR(phi(i, j, k), 0.0, 1e-4);
}

TEST(Dense3DSolverTest, StressTestModerateGrid)
{
    const size_t n = 15;
    const double dx = 0.5;
    Dense3DSolver solver((n - 1) * dx, (n - 1) * dx, (n - 1) * dx, dx);
    ASSERT_EQ(solver.nx(), n);

    auto rhs = make_zero_rhs(solver.nx(), solver.ny(), solver.nz());
    solver.setRHS(rhs);

    std::array<double, 6> bc = {25.0, 40.0, 55.0, 70.0, 85.0, 100.0};
    solver.setBoundary(bc);

    solver.solve(4000, 1e-5);

    const auto &phi = solver.solution();
    double center = phi(n / 2, n / 2, n / 2);
    EXPECT_GT(center, 25.0);
    EXPECT_LT(center, 100.0);
}


// ============================================================
// Golden Test 3D: T = sin(pi x) sin(pi y) sinh(pi z)/sinh(pi L)
// ============================================================
TEST(GoldenTest, Analytic3D)
{
    const size_t NX = 20, NY = 20, NZ = 20;
    const double L = 1.0;

    double dx = L / (NX - 1);
    double dy = L / (NY - 1);
    double dz = L / (NZ - 1);

    Dense3DSolver solver(L, L, L, dx);
    ASSERT_EQ(solver.nx(), NX);
    auto rhs = make_zero_rhs(solver.nx(), solver.ny(), solver.nz());
    solver.setRHS(rhs);

    std::array<double, 6> bc = {0, 0, 0, 0, 0, 0};
    solver.setBoundary(bc);

    solver.solve(10000, 1e-6);

    const auto &phi = solver.solution();
    for (size_t k = 0; k < NZ; k++)
        for (size_t j = 0; j < NY; j++)
            for (size_t i = 0; i < NX; i++)
            {
                double xv = i * dx, yv = j * dy, zv = k * dz;
                double expected = sin(M_PI * xv) *
                                  sin(M_PI * yv) *
                                  sinh(M_PI * zv) / sinh(M_PI * L);

                EXPECT_NEAR(phi(i, j, k), expected, 1.0);
            }
}

//
// ============================================================
// Test main
// ============================================================
int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
