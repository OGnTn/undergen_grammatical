#ifndef GRID_PARALLEL_H
#define GRID_PARALLEL_H

#include <algorithm>
#include <cstdint>
#include <thread>
#include <vector>

namespace godot {

inline int grid_parallel_worker_count(int z_count, int64_t total_cells) {
    if (z_count <= 0 || total_cells < 262144) {
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

    // Thread startup costs dominate small grids; keep those single-threaded.
    if (worker_count == 1 || total_cells < 262144) {
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
        worker.join();
    }
}

} // namespace godot

#endif // GRID_PARALLEL_H
