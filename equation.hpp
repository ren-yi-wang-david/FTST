#pragma once
#include <memory>
#include <cstddef>
#include <algorithm>
#include <stdexcept>
#include"memory_tracker.hpp"
template <typename T, typename Alloc = TrackingAllocator<T>>
class Matrix3D

{
public:
    using value_type = T;
    using allocator_type = Alloc;
    using traits = std::allocator_traits<allocator_type>;

    // -------------------------
    // Constructor
    // -------------------------
    Matrix3D(size_t nx, size_t ny, size_t nz,
             const allocator_type &alloc = allocator_type())
        : nx_(nx), ny_(ny), nz_(nz), alloc_(alloc)
    {
        allocate_and_construct();
    }

    // -------------------------
    // Destructor (Rule of Five #1)
    // -------------------------
    ~Matrix3D()
    {
        destroy_and_deallocate();
    }

    // -------------------------
    // Copy Constructor (Rule of Five #2)
    // -------------------------
    Matrix3D(const Matrix3D &other)
        : nx_(other.nx_), ny_(other.ny_), nz_(other.nz_),
          alloc_(traits::select_on_container_copy_construction(other.alloc_))
    {
        allocate_and_construct();
        std::uninitialized_copy(other.data_, other.data_ + total_size(), data_);
    }

    // -------------------------
    // Copy Assignment (Rule of Five #3)
    // -------------------------
    Matrix3D &operator=(const Matrix3D &other)
    {
        if (this == &other)
            return *this;

        // Strong exception guarantee:
        // copy first, then swap
        Matrix3D tmp(other);
        swap(tmp);
        return *this;
    }

    // -------------------------
    // Move Constructor (Rule of Five #4)
    // -------------------------
    Matrix3D(Matrix3D &&other) noexcept
        : nx_(other.nx_), ny_(other.ny_), nz_(other.nz_),
          alloc_(std::move(other.alloc_)), data_(other.data_)
    {
        other.data_ = nullptr;
        other.nx_ = other.ny_ = other.nz_ = 0;
    }

    // -------------------------
    // Move Assignment (Rule of Five #5)
    // -------------------------
    Matrix3D &operator=(Matrix3D &&other) noexcept
    {
        if (this == &other)
            return *this;

        destroy_and_deallocate();

        nx_ = other.nx_;
        ny_ = other.ny_;
        nz_ = other.nz_;
        data_ = other.data_;
        alloc_ = std::move(other.alloc_);

        other.data_ = nullptr;
        other.nx_ = other.ny_ = other.nz_ = 0;

        return *this;
    }

    // -------------------------
    // Element Access
    // -------------------------
    T &operator()(size_t i, size_t j, size_t k) noexcept
    {
        return data_[index(i, j, k)];
    }

    const T &operator()(size_t i, size_t j, size_t k) const noexcept
    {
        return data_[index(i, j, k)];
    }

    // -------------------------
    // Basic Info
    // -------------------------
    size_t nx() const noexcept { return nx_; }
    size_t ny() const noexcept { return ny_; }
    size_t nz() const noexcept { return nz_; }
    size_t total_size() const noexcept { return nx_ * ny_ * nz_; }

    T *data() noexcept { return data_; }
    const T *data() const noexcept { return data_; }

    // -------------------------
    // Swap (Strong Guarantee Tool)
    // -------------------------
    void swap(Matrix3D &other) noexcept
    {
        std::swap(nx_, other.nx_);
        std::swap(ny_, other.ny_);
        std::swap(nz_, other.nz_);
        std::swap(data_, other.data_);
        std::swap(alloc_, other.alloc_);
    }

private:
    size_t index(size_t i, size_t j, size_t k) const noexcept
    {
        return (i * ny_ + j) * nz_ + k;
    }

    // -------------------------
    // Memory Ops
    // -------------------------
    void allocate_and_construct()
    {
        if (total_size() == 0)
        {
            data_ = nullptr;
            return;
        }

        data_ = traits::allocate(alloc_, total_size());
        for (size_t i = 0; i < total_size(); ++i)
            traits::construct(alloc_, data_ + i, T{});
    }

    void destroy_and_deallocate()
    {
        if (!data_)
            return;

        for (size_t i = 0; i < total_size(); ++i)
            traits::destroy(alloc_, data_ + i);

        traits::deallocate(alloc_, data_, total_size());
        data_ = nullptr;
    }

private:
    size_t nx_, ny_, nz_;
    allocator_type alloc_;
    T *data_;
};
