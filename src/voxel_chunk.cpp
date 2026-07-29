#include "voxel_chunk.h"

namespace godot {

VoxelChunk::VoxelChunk() : chunk_coord(0, 0, 0), state(CHUNK_STATE_UNIFORM_AIR) {}

VoxelChunk::VoxelChunk(const Vector3i &coord, VoxelChunkState initial_state, float air_val)
    : chunk_coord(coord), state(initial_state), default_air_value(air_val) {}

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
    density_data.resize(CHUNK_CELL_COUNT, fill_density);
    material_data.resize(CHUNK_CELL_COUNT, fill_material);
    zone_data.resize(CHUNK_CELL_COUNT, fill_zone);
    is_dirty = true;
}

} // namespace godot
