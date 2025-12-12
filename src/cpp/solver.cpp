#include <ftst/solver.hpp>
#include <cmath>
#include <array>
#include <immintrin.h>
#include <thread>
#include <vector>
#include <algorithm>

// ============================================================
// Small helper
// ============================================================
namespace
{
    size_t to_steps(double length, double spacing)
    {
        return static_cast<size_t>(std::llround(length / spacing)) + 1;
    }
}

// ============================================================
// Constructor
// ============================================================
Dense3DSolver::Dense3DSolver(size_t nx, size_t ny, size_t nz, double dx)
    : nx_(nx), ny_(ny), nz_(nz), dx_(dx),
      phi_(nx, ny, nz),
      rhs_(nx, ny, nz)
{
    for (size_t i = 0; i < nx_; ++i)
        for (size_t j = 0; j < ny_; ++j)
            for (size_t k = 0; k < nz_; ++k)
                phi_(i, j, k) = 0.0;

    // Best SOR relaxation weight
    double rx = std::cos(M_PI / nx_);
    double ry = std::cos(M_PI / ny_);
    double rz = std::cos(M_PI / nz_);
    double rho = (rx + ry + rz) / 3.0;
    omega = 2.0 / (1.0 + std::sqrt(1.0 - rho * rho));

    flattenMemory();
    buildColoring();
}

Dense3DSolver::Dense3DSolver(double lx, double ly, double lz, double dx)
    : Dense3DSolver(to_steps(lx, dx), to_steps(ly, dx), to_steps(lz, dx), dx)
{
}

// ============================================================
// Boundary check (O(1))
// ============================================================
inline bool Dense3DSolver::isBoundary(size_t idx) const
{
    size_t k = idx % nz_;
    size_t j = (idx / nz_) % ny_;
    size_t i = idx / (ny_ * nz_);

    return (i == 0 || i == nx_ - 1 ||
            j == 0 || j == ny_ - 1 ||
            k == 0 || k == nz_ - 1);
}

// ============================================================
// Set RHS
// ============================================================
void Dense3DSolver::setRHS(const Matrix3DDouble &rhs)
{
    for (size_t i = 0; i < nx_; ++i)
        for (size_t j = 0; j < ny_; ++j)
            for (size_t k = 0; k < nz_; ++k)
            {
                rhs_(i, j, k) = rhs(i, j, k);
                rhs1D_[i * ny_ * nz_ + j * nz_ + k] = rhs(i, j, k);
            }
}

// ============================================================
// Set boundary
// ============================================================
void Dense3DSolver::setBoundary(const std::array<double, 6> &bc)
{
    // x-min and x-max
    for (size_t j = 0; j < ny_; ++j)
        for (size_t k = 0; k < nz_; ++k)
        {
            phi_(0, j, k) = bc[0];
            phi_(nx_ - 1, j, k) = bc[1];

            phi1D_[0 * ny_ * nz_ + j * nz_ + k] = bc[0];
            phi1D_[(nx_ - 1) * ny_ * nz_ + j * nz_ + k] = bc[1];
        }

    // y-min and y-max
    for (size_t i = 0; i < nx_; ++i)
        for (size_t k = 0; k < nz_; ++k)
        {
            phi_(i, 0, k) = bc[2];
            phi_(i, ny_ - 1, k) = bc[3];

            phi1D_[i * ny_ * nz_ + 0 * nz_ + k] = bc[2];
            phi1D_[i * ny_ * nz_ + (ny_ - 1) * nz_ + k] = bc[3];
        }

    // z-min and z-max
    for (size_t i = 0; i < nx_; ++i)
        for (size_t j = 0; j < ny_; ++j)
        {
            phi_(i, j, 0) = bc[4];
            phi_(i, j, nz_ - 1) = bc[5];

            phi1D_[i * ny_ * nz_ + j * nz_ + 0] = bc[4];
            phi1D_[i * ny_ * nz_ + j * nz_ + (nz_ - 1)] = bc[5];
        }
}

// ============================================================
// Flatten 3D into 1D
// ============================================================
void Dense3DSolver::flattenMemory()
{
    size_t N = nx_ * ny_ * nz_;
    phi1D_.assign(N, 0.0);
    rhs1D_.assign(N, 0.0);
}

// ============================================================
// Build red/black lists and skip boundaries
// ============================================================
void Dense3DSolver::buildColoring()
{
    red_list_.clear();
    black_list_.clear();

    const size_t Ti = 8, Tj = 8, Tk = 8;

    for (size_t ii = 1; ii < nx_ - 1; ii += Ti)
        for (size_t jj = 1; jj < ny_ - 1; jj += Tj)
            for (size_t kk = 1; kk < nz_ - 1; kk += Tk)
            {
                size_t ie = std::min(ii + Ti, nx_ - 1);
                size_t je = std::min(jj + Tj, ny_ - 1);
                size_t ke = std::min(kk + Tk, nz_ - 1);

                for (size_t i = ii; i < ie; ++i)
                    for (size_t j = jj; j < je; ++j)
                        for (size_t k = kk; k < ke; ++k)
                        {
                            size_t idx = i * ny_ * nz_ + j * nz_ + k;

                            if (isBoundary(idx))
                                continue;

                            if (((i + j + k) & 1) == 0)
                                red_list_.push_back(idx);
                            else
                                black_list_.push_back(idx);
                        }
            }
}

// ============================================================
// Solve using Red-Black SOR + Threads + AVX2
// ============================================================
void Dense3DSolver::solve(size_t maxIter, double tol)
{
    const double dx2 = dx_ * dx_;
    const double w = omega;

    const size_t ox = ny_ * nz_;
    const size_t oy = nz_;
    const size_t oz = 1;

    const __m256d vdx2 = _mm256_set1_pd(dx2);
    const __m256d vw = _mm256_set1_pd(w);
    const __m256d v6 = _mm256_set1_pd(6.0);

    const __m256i vox = _mm256_set1_epi64x(ox);
    const __m256i voy = _mm256_set1_epi64x(oy);
    const __m256i voz = _mm256_set1_epi64x(oz);

    int numThreads = std::thread::hardware_concurrency();
    if (numThreads <= 0)
        numThreads = 8;

    auto process_segment = [&](const std::vector<size_t> &list,
                               size_t start, size_t end,
                               double &localMax)
    {
        double myMax = 0.0;

        size_t n = end - start;
        size_t vecEnd = start + (n / 4) * 4;

        // SIMD section
        for (size_t t = start; t < vecEnd; t += 4)
        {
            __m256i vidx = _mm256_set_epi64x(
                list[t + 3], list[t + 2], list[t + 1], list[t]);

            __m256d vphi = _mm256_i64gather_pd(phi1D_.data(), vidx, 8);

            __m256d vpx = _mm256_i64gather_pd(phi1D_.data(), _mm256_add_epi64(vidx, vox), 8);
            __m256d vmx = _mm256_i64gather_pd(phi1D_.data(), _mm256_sub_epi64(vidx, vox), 8);

            __m256d vpy = _mm256_i64gather_pd(phi1D_.data(), _mm256_add_epi64(vidx, voy), 8);
            __m256d vmy = _mm256_i64gather_pd(phi1D_.data(), _mm256_sub_epi64(vidx, voy), 8);

            __m256d vpz = _mm256_i64gather_pd(phi1D_.data(), _mm256_add_epi64(vidx, voz), 8);
            __m256d vmz = _mm256_i64gather_pd(phi1D_.data(), _mm256_sub_epi64(vidx, voz), 8);

            __m256d vrhs = _mm256_i64gather_pd(rhs1D_.data(), vidx, 8);

            __m256d vsum = _mm256_add_pd(
                _mm256_add_pd(_mm256_add_pd(vpx, vmx), _mm256_add_pd(vpy, vmy)),
                _mm256_add_pd(vpz, vmz));

            __m256d vgs = _mm256_div_pd(
                _mm256_sub_pd(vsum, _mm256_mul_pd(vdx2, vrhs)), v6);

            __m256d vdiff = _mm256_sub_pd(vgs, vphi);
            __m256d vnew = _mm256_add_pd(vphi, _mm256_mul_pd(vw, vdiff));

            alignas(32) double buf[4];
            _mm256_store_pd(buf, vnew);

            size_t id0 = list[t];
            size_t id1 = list[t + 1];
            size_t id2 = list[t + 2];
            size_t id3 = list[t + 3];

            if (!isBoundary(id0))
                phi1D_[id0] = buf[0];
            if (!isBoundary(id1))
                phi1D_[id1] = buf[1];
            if (!isBoundary(id2))
                phi1D_[id2] = buf[2];
            if (!isBoundary(id3))
                phi1D_[id3] = buf[3];

            alignas(32) double dBuf[4];
            _mm256_store_pd(dBuf, vdiff);
            for (double x : dBuf)
                myMax = std::max(myMax, std::abs(x));
        }

        // scalar tail
        for (size_t t = vecEnd; t < end; ++t)
        {
            size_t idx = list[t];
            double old = phi1D_[idx];

            double gs =
                (phi1D_[idx + ox] + phi1D_[idx - ox] +
                 phi1D_[idx + oy] + phi1D_[idx - oy] +
                 phi1D_[idx + oz] + phi1D_[idx - oz] -
                 dx2 * rhs1D_[idx]) /
                6.0;

            double newv = old + w * (gs - old);

            if (!isBoundary(idx))
                phi1D_[idx] = newv;

            myMax = std::max(myMax, std::abs(newv - old));
        }

        localMax = myMax;
    };

    // ============================
    // Main SOR iteration
    // ============================
    for (size_t iter = 0; iter < maxIter; ++iter)
    {
        double maxDiff = 0.0;

        // ------- RED -------
        {
            size_t N = red_list_.size();
            size_t chunk = (N + numThreads - 1) / numThreads;

            std::vector<std::thread> threads;
            std::vector<double> localMax(numThreads, 0.0);

            for (int t = 0; t < numThreads; ++t)
            {
                size_t start = t * chunk;
                size_t end = std::min(start + chunk, N);
                if (start >= end)
                    break;

                threads.emplace_back(process_segment,
                                     std::ref(red_list_),
                                     start, end,
                                     std::ref(localMax[t]));
            }
            for (auto &th : threads)
                th.join();
            for (double v : localMax)
                maxDiff = std::max(maxDiff, v);
        }

        // ------- BLACK -------
        {
            size_t N = black_list_.size();
            size_t chunk = (N + numThreads - 1) / numThreads;

            std::vector<std::thread> threads;
            std::vector<double> localMax(numThreads, 0.0);

            for (int t = 0; t < numThreads; ++t)
            {
                size_t start = t * chunk;
                size_t end = std::min(start + chunk, N);
                if (start >= end)
                    break;

                threads.emplace_back(process_segment,
                                     std::ref(black_list_),
                                     start, end,
                                     std::ref(localMax[t]));
            }
            for (auto &th : threads)
                th.join();
            for (double v : localMax)
                maxDiff = std::max(maxDiff, v);
        }

        if (maxDiff < tol)
            break;
    }

    // copy back into 3D
    for (size_t i = 0; i < nx_; ++i)
        for (size_t j = 0; j < ny_; ++j)
            for (size_t k = 0; k < nz_; ++k)
                phi_(i, j, k) = phi1D_[i * ny_ * nz_ + j * nz_ + k];
}
