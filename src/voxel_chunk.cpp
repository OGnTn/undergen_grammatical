#include "voxel_chunk.h"

namespace godot {

VoxelChunk::VoxelChunk() : chunk_coord(0, 0, 0), state(CHUNK_STATE_UNIFORM_AIR), octree(CHUNK_SIZE, 1.0f) {}

VoxelChunk::VoxelChunk(const Vector3i &coord, VoxelChunkState initial_state, float air_val)
    : chunk_coord(coord), state(initial_state), default_air_value(air_val),
      octree(CHUNK_SIZE, initial_state == CHUNK_STATE_UNIFORM_SOLID ? -1.0f : air_val) {}

VoxelChunk::~VoxelChunk() {
    density_data.clear();
    material_data.clear();
    zone_data.clear();
}

void VoxelChunk::allocate_data(float fill_density, uint8_t fill_material, int32_t fill_zone) {
    if (state == CHUNK_STATE_ALLOCATED) {
        return;
    }
    state = CHUNK_STATE_ALLOCATED;
    octree.initialize(CHUNK_SIZE, fill_density, fill_material, fill_zone);
    is_dirty = true;
}

void VoxelChunk::ensure_direct_ptr_allocated() {
    if (!direct_ptr_allocated) {
        if (state != CHUNK_STATE_ALLOCATED) {
            allocate_data(state == CHUNK_STATE_UNIFORM_SOLID ? default_solid_value : default_air_value);
        }
        density_data.resize(CHUNK_CELL_COUNT);
        material_data.resize(CHUNK_CELL_COUNT);
        zone_data.resize(CHUNK_CELL_COUNT);

        for (int z = 0; z < CHUNK_SIZE; ++z) {
            for (int y = 0; y < CHUNK_SIZE; ++y) {
                for (int x = 0; x < CHUNK_SIZE; ++x) {
                    int idx = get_local_index(x, y, z);
                    Vector3i p(x, y, z);
                    density_data[idx] = octree.get_density(p);
                    material_data[idx] = octree.get_material(p);
                    zone_data[idx] = octree.get_zone(p);
                }
            }
        }
        direct_ptr_allocated = true;
    }
}

} // namespace godot
