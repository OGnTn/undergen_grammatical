#ifndef GRID_PARALLEL_H
#define GRID_PARALLEL_H

#include <algorithm>
#include <cstdint>
#include <thread>
#include <vector>
#include <functional>

namespace godot {

inline int grid_parallel_worker_count(int z_count, int64_t total_cells) {
    if (z_count <= 0) {
        return 1;
    }
    // Allow parallel execution for grids with > 32k cells
    if (total_cells < 32768) {
        return 1;
    }

    unsigned int hw_threads = std::thread::hardware_concurrency();
    int worker_count = (int)(hw_threads > 0 ? hw_threads : 4);
    return std::max(1, std::min(worker_count, z_count));
}

template <typename Func>
void parallel_for_z(int z_count, int64_t total_cells, Func func) {
    if (z_count <= 0) {
        return;
    }

    int worker_count = grid_parallel_worker_count(z_count, total_cells);

    if (worker_count == 1) {
        func(0, 0, z_count);
        return;
    }

    std::vector<std::thread> workers;
    workers.reserve(worker_count);

    int z_begin = 0;
    for (int worker = 0; worker < worker_count; ++worker) {
        int remaining_z = z_count - z_begin;
        int remaining_workers = worker_count - worker;
        int slice_count = (remaining_z + remaining_workers - 1) / remaining_workers;
        int z_end = z_begin + slice_count;

        workers.emplace_back([=, &func]() {
            func(worker, z_begin, z_end);
        });

        z_begin = z_end;
    }

    for (std::thread &worker : workers) {
        if (worker.joinable()) {
            worker.join();
        }
    }
}

template <typename ChunkType, typename Func>
void parallel_for_chunks(const std::vector<ChunkType*> &chunks, Func func) {
    size_t count = chunks.size();
    if (count == 0) return;

    unsigned int hw_threads = std::thread::hardware_concurrency();
    int worker_count = std::max(1, std::min((int)hw_threads, (int)count));

    if (worker_count == 1) {
        func(0, 0, (int)count);
        return;
    }

    std::vector<std::thread> workers;
    workers.reserve(worker_count);

    int begin = 0;
    for (int worker = 0; worker < worker_count; ++worker) {
        int remaining = (int)count - begin;
        int remaining_workers = worker_count - worker;
        int chunk_slice = (remaining + remaining_workers - 1) / remaining_workers;
        int end = begin + chunk_slice;

        workers.emplace_back([=, &func]() {
            func(worker, begin, end);
        });

        begin = end;
    }

    for (std::thread &worker : workers) {
        if (worker.joinable()) {
            worker.join();
        }
    }
}

} // namespace godot

#endif // GRID_PARALLEL_H
