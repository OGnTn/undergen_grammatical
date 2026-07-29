#ifndef VOXEL_CHUNK_H
#define VOXEL_CHUNK_H

#include <vector>
#include <cstdint>
#include <godot_cpp/variant/vector3i.hpp>

namespace godot {

constexpr int CHUNK_SIZE = 16;
constexpr int CHUNK_SIZE_SHIFT = 4;
constexpr int CHUNK_SIZE_MASK = 15;
constexpr int CHUNK_CELL_COUNT = CHUNK_SIZE * CHUNK_SIZE * CHUNK_SIZE; // 4096

enum VoxelChunkState {
    CHUNK_STATE_UNIFORM_AIR = 0,
    CHUNK_STATE_UNIFORM_SOLID = 1,
    CHUNK_STATE_ALLOCATED = 2
};

class VoxelChunk {
private:
    Vector3i chunk_coord;
    VoxelChunkState state = CHUNK_STATE_UNIFORM_AIR;
    float default_air_value = 1.0f;
    float default_solid_value = -1.0f;

    std::vector<float> density_data;
    std::vector<uint8_t> material_data;
    std::vector<int32_t> zone_data;

    bool is_dirty = false;

public:
    VoxelChunk();
    VoxelChunk(const Vector3i &coord, VoxelChunkState initial_state = CHUNK_STATE_UNIFORM_AIR, float air_val = 1.0f);
    ~VoxelChunk();

    Vector3i get_chunk_coord() const { return chunk_coord; }
    VoxelChunkState get_state() const { return state; }
    bool is_allocated() const { return state == CHUNK_STATE_ALLOCATED; }
    bool get_is_dirty() const { return is_dirty; }
    void set_is_dirty(bool dirty) { is_dirty = dirty; }

    void allocate_data(float fill_density = 1.0f, uint8_t fill_material = 0, int32_t fill_zone = 0);

    static inline int get_local_index(int x, int y, int z) {
        return x + CHUNK_SIZE * (y + CHUNK_SIZE * z);
    }

    inline float get_density(int x, int y, int z) const {
        if (state == CHUNK_STATE_ALLOCATED) {
            return density_data[get_local_index(x, y, z)];
        }
        return (state == CHUNK_STATE_UNIFORM_SOLID) ? default_solid_value : default_air_value;
    }

    inline void set_density(int x, int y, int z, float val) {
        if (state != CHUNK_STATE_ALLOCATED) {
            allocate_data(state == CHUNK_STATE_UNIFORM_SOLID ? default_solid_value : default_air_value);
        }
        density_data[get_local_index(x, y, z)] = val;
        is_dirty = true;
    }

    inline uint8_t get_material(int x, int y, int z) const {
        if (state == CHUNK_STATE_ALLOCATED) {
            return material_data[get_local_index(x, y, z)];
        }
        return 0;
    }

    inline void set_material(int x, int y, int z, uint8_t mat_id) {
        if (state != CHUNK_STATE_ALLOCATED) {
            allocate_data(state == CHUNK_STATE_UNIFORM_SOLID ? default_solid_value : default_air_value);
        }
        material_data[get_local_index(x, y, z)] = mat_id;
        is_dirty = true;
    }

    inline int32_t get_zone(int x, int y, int z) const {
        if (state == CHUNK_STATE_ALLOCATED) {
            return zone_data[get_local_index(x, y, z)];
        }
        return 0;
    }

    inline void set_zone(int x, int y, int z, int32_t zone_id) {
        if (state != CHUNK_STATE_ALLOCATED) {
            allocate_data(state == CHUNK_STATE_UNIFORM_SOLID ? default_solid_value : default_air_value);
        }
        zone_data[get_local_index(x, y, z)] = zone_id;
        is_dirty = true;
    }

    // Direct pointer access for fast C++ operations when allocated
    float* get_density_ptr() { return density_data.empty() ? nullptr : density_data.data(); }
    const float* get_density_ptr() const { return density_data.empty() ? nullptr : density_data.data(); }

    uint8_t* get_material_ptr() { return material_data.empty() ? nullptr : material_data.data(); }
    const uint8_t* get_material_ptr() const { return material_data.empty() ? nullptr : material_data.data(); }

    int32_t* get_zone_ptr() { return zone_data.empty() ? nullptr : zone_data.data(); }
    const int32_t* get_zone_ptr() const { return zone_data.empty() ? nullptr : zone_data.data(); }
};

} // namespace godot

#endif // VOXEL_CHUNK_H
