#pragma once
#include <memory>
#include <cstddef>
#include <atomic>
#include <iostream>
#include <algorithm>

struct MemoryTracker
{
    inline static std::atomic<size_t> bytes_allocated{0};
    inline static std::atomic<size_t> bytes_freed{0};
    inline static std::atomic<size_t> peak{0};

    static void report()
    {
        std::cout << "== == == == == == == == == == == == == ==" << std::endl;
        std::cout << "Allocated: " << bytes_allocated.load() << " bytes\n";
        std::cout << "Freed:     " << bytes_freed.load() << " bytes\n";
        std::cout << "Peak:      " << peak.load() << " bytes\n";
        std::cout << "== == == == == == == == == == == == == ==" << std::endl;
    }
};

template <typename T>
class TrackingAllocator
{
public:
    using value_type = T;

    TrackingAllocator() noexcept {}
    template <class U>
    TrackingAllocator(const TrackingAllocator<U> &) noexcept {}

    T *allocate(std::size_t n)
    {
        std::size_t bytes = n * sizeof(T);
        MemoryTracker::bytes_allocated += bytes;

        size_t current = MemoryTracker::bytes_allocated - MemoryTracker::bytes_freed;
        MemoryTracker::peak = std::max(MemoryTracker::peak.load(), current);

        return static_cast<T *>(::operator new(bytes));
    }

    void deallocate(T *p, std::size_t n) noexcept
    {
        std::size_t bytes = n * sizeof(T);
        MemoryTracker::bytes_freed += bytes;
        ::operator delete(p);
    }
};

template <class T, class U>
bool operator==(const TrackingAllocator<T> &, const TrackingAllocator<U> &) { return true; }

template <class T, class U>
bool operator!=(const TrackingAllocator<T> &, const TrackingAllocator<U> &) { return false; }
