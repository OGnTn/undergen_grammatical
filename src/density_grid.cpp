// src/density_grid.cpp
#include "density_grid.h"

#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/variant/utility_functions.hpp>

namespace godot {

DensityGrid::DensityGrid() {}

DensityGrid::~DensityGrid() {
    clear_chunks();
}

void DensityGrid::clear_chunks() {
    for (auto &pair : chunks) {
        if (pair.second) {
            delete pair.second;
        }
    }
    chunks.clear();
    clear_hermite_edges();
}

Vector3i DensityGrid::world_to_chunk_coord(const Vector3i &pos) {
    return Vector3i(pos.x >> CHUNK_SIZE_SHIFT, pos.y >> CHUNK_SIZE_SHIFT, pos.z >> CHUNK_SIZE_SHIFT);
}

Vector3i DensityGrid::world_to_local_coord(const Vector3i &pos) {
    return Vector3i(pos.x & CHUNK_SIZE_MASK, pos.y & CHUNK_SIZE_MASK, pos.z & CHUNK_SIZE_MASK);
}

VoxelChunk* DensityGrid::get_chunk(const Vector3i &chunk_coord, bool create_if_missing) {
    auto it = chunks.find(chunk_coord);
    if (it != chunks.end()) {
        return it->second;
    }
    if (!create_if_missing) {
        return nullptr;
    }

    VoxelChunkState initial_state = (default_initial_value <= -0.5f) ? CHUNK_STATE_UNIFORM_SOLID : CHUNK_STATE_UNIFORM_AIR;
    VoxelChunk *chunk = new VoxelChunk(chunk_coord, initial_state, default_initial_value);
    chunks[chunk_coord] = chunk;
    return chunk;
}

VoxelChunk* DensityGrid::get_chunk_for_voxel(const Vector3i &pos, bool create_if_missing) {
    Vector3i c_coord = world_to_chunk_coord(pos);
    return get_chunk(c_coord, create_if_missing);
}

void DensityGrid::_bind_methods() {
    ClassDB::bind_method(D_METHOD("initialize_grid", "size_x", "size_y", "size_z", "default_value", "default_material_index"), &DensityGrid::initialize_grid, DEFVAL(1.0), DEFVAL(0));
    ClassDB::bind_method(D_METHOD("is_valid_position", "pos"), &DensityGrid::is_valid_position);
    ClassDB::bind_method(D_METHOD("get_index", "pos"), &DensityGrid::get_index);
    ClassDB::bind_method(D_METHOD("set_cell", "pos", "value"), &DensityGrid::set_cell);
    ClassDB::bind_method(D_METHOD("get_cell", "pos", "default_value"), &DensityGrid::get_cell, DEFVAL(1.0));
    ClassDB::bind_method(D_METHOD("get_grid_dimensions"), &DensityGrid::get_grid_dimensions);

    ClassDB::bind_method(D_METHOD("set_world_density_grid", "grid"), &DensityGrid::set_world_density_grid);
    ClassDB::bind_method(D_METHOD("get_world_density_grid"), &DensityGrid::get_world_density_grid);
    ADD_PROPERTY(PropertyInfo(Variant::PACKED_FLOAT32_ARRAY, "world_density_grid"), "set_world_density_grid", "get_world_density_grid");

    ClassDB::bind_method(D_METHOD("set_grid_size_x", "x"), &DensityGrid::set_grid_size_x);
    ClassDB::bind_method(D_METHOD("get_grid_size_x"), &DensityGrid::get_grid_size_x);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "grid_size_x"), "set_grid_size_x", "get_grid_size_x");

    ClassDB::bind_method(D_METHOD("set_grid_size_y", "y"), &DensityGrid::set_grid_size_y);
    ClassDB::bind_method(D_METHOD("get_grid_size_y"), &DensityGrid::get_grid_size_y);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "grid_size_y"), "set_grid_size_y", "get_grid_size_y");

    ClassDB::bind_method(D_METHOD("set_grid_size_z", "z"), &DensityGrid::set_grid_size_z);
    ClassDB::bind_method(D_METHOD("get_grid_size_z"), &DensityGrid::get_grid_size_z);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "grid_size_z"), "set_grid_size_z", "get_grid_size_z");

    ClassDB::bind_method(D_METHOD("set_surface_threshold", "threshold"), &DensityGrid::set_surface_threshold);
    ClassDB::bind_method(D_METHOD("get_surface_threshold"), &DensityGrid::get_surface_threshold);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "surface_threshold"), "set_surface_threshold", "get_surface_threshold");

    ClassDB::bind_method(D_METHOD("set_material_id", "pos", "material_index"), &DensityGrid::set_material_id);
    ClassDB::bind_method(D_METHOD("get_material_id", "pos"), &DensityGrid::get_material_id);

    ClassDB::bind_method(D_METHOD("set_world_material_grid", "grid"), &DensityGrid::set_world_material_grid);
    ClassDB::bind_method(D_METHOD("get_world_material_grid"), &DensityGrid::get_world_material_grid);
    ADD_PROPERTY(PropertyInfo(Variant::PACKED_BYTE_ARRAY, "world_material_grid"), "set_world_material_grid", "get_world_material_grid");

    ClassDB::bind_method(D_METHOD("get_zone_name_by_id", "id"), &DensityGrid::get_zone_name_by_id);
    ClassDB::bind_method(D_METHOD("get_zone_at", "pos"), &DensityGrid::get_zone_at);
}

void DensityGrid::initialize_grid(int size_x, int size_y, int size_z, float default_value, int default_material_index) {
    if (size_x <= 0 || size_y <= 0 || size_z <= 0) {
        UtilityFunctions::printerr("DensityGrid.initialize_grid: Invalid dimensions provided (", size_x, ", ", size_y, ", ", size_z, ").");
        grid_size_x = 0;
        grid_size_y = 0;
        grid_size_z = 0;
        world_density_grid.clear();
        world_material_grid.clear();
        world_zone_grid.clear();
        clear_chunks();
        return;
    }

    grid_size_x = size_x;
    grid_size_y = size_y;
    grid_size_z = size_z;
    default_initial_value = default_value;
    default_material_idx = default_material_index;

    int64_t total_grid_points = (int64_t)grid_size_x * grid_size_y * grid_size_z;
    world_density_grid.resize(total_grid_points);
    world_density_grid.fill(default_value);

    world_material_grid.resize(total_grid_points);
    world_material_grid.fill(default_material_index);

    world_zone_grid.resize(total_grid_points);
    world_zone_grid.fill(0);

    clear_chunks();

    zone_id_map.clear();
    zone_id_map.push_back("undefined");

    UtilityFunctions::print("Initialized density grid of size ", grid_size_x, " x ", grid_size_y, " x ", grid_size_z, " = ", total_grid_points, " points with default value ", default_value);
}

bool DensityGrid::is_valid_position(const Vector3i &pos) const {
    return pos.x >= 0 && pos.x < grid_size_x &&
           pos.y >= 0 && pos.y < grid_size_y &&
           pos.z >= 0 && pos.z < grid_size_z;
}

int DensityGrid::get_index(const Vector3i &pos) const {
    if (!is_valid_position(pos)) {
        return -1;
    }
    int64_t index = (int64_t)pos.x + (int64_t)grid_size_x * (pos.y + (int64_t)grid_size_y * pos.z);
    if (index < 0 || index >= world_density_grid.size()) {
        return -1;
    }
    return (int)index;
}

bool DensityGrid::set_cell(const Vector3i &pos, float value) {
    int index = get_index(pos);
    if (index != -1 && index < world_density_grid.size()) {
        world_density_grid[index] = value;

        Vector3i l_pos = world_to_local_coord(pos);
        VoxelChunk *chunk = get_chunk_for_voxel(pos, true);
        if (chunk) {
            chunk->set_density(l_pos.x, l_pos.y, l_pos.z, value);
        }
        return true;
    }
    return false;
}

float DensityGrid::get_cell(const Vector3i &pos, float default_value) const {
    if (grid_size_x <= 0 || grid_size_y <= 0 || grid_size_z <= 0) {
        return default_value;
    }

    if (pos.y >= grid_size_y) {
        return -1.0f; // Sky above the grid is empty/air
    }
    if (pos.y < 0) {
        return 1.0f; // Bedrock below the grid is solid
    }

    int clamped_x = pos.x < 0 ? 0 : (pos.x >= grid_size_x ? grid_size_x - 1 : pos.x);
    int clamped_z = pos.z < 0 ? 0 : (pos.z >= grid_size_z ? grid_size_z - 1 : pos.z);

    Vector3i clamped_pos(clamped_x, pos.y, clamped_z);
    int index = get_index(clamped_pos);
    if (index >= 0 && index < world_density_grid.size()) {
        return world_density_grid[index];
    }
    return default_value;
}

Vector3i DensityGrid::get_grid_dimensions() const {
    return Vector3i(grid_size_x, grid_size_y, grid_size_z);
}

void DensityGrid::sync_flat_cache_if_dirty() const {
    // Kept for backward compatibility
}

void DensityGrid::set_world_density_grid(const PackedFloat32Array &p_grid) {
    world_density_grid = p_grid;
}

PackedFloat32Array DensityGrid::get_world_density_grid() const {
    return world_density_grid;
}

PackedFloat32Array &DensityGrid::get_density_data_rw() {
    return world_density_grid;
}

const PackedFloat32Array &DensityGrid::get_density_data() const {
    return world_density_grid;
}

PackedByteArray &DensityGrid::get_material_data_rw() {
    return world_material_grid;
}

const PackedByteArray &DensityGrid::get_material_data() const {
    return world_material_grid;
}

PackedInt32Array &DensityGrid::get_zone_data_rw() {
    return world_zone_grid;
}

const PackedInt32Array &DensityGrid::get_zone_data() const {
    return world_zone_grid;
}

int64_t DensityGrid::get_total_cell_count() const {
    return (int64_t)grid_size_x * grid_size_y * grid_size_z;
}

int DensityGrid::get_index_unchecked(int x, int y, int z) const {
    return x + grid_size_x * (y + grid_size_y * z);
}

void DensityGrid::set_material_id(const Vector3i &pos, int material_index) {
    int index = get_index(pos);
    if (index != -1 && index < world_material_grid.size()) {
        if (material_index < 0) material_index = 0;
        if (material_index > 255) material_index = 255;
        world_material_grid[index] = (uint8_t)material_index;

        Vector3i l_pos = world_to_local_coord(pos);
        VoxelChunk *chunk = get_chunk_for_voxel(pos, true);
        if (chunk) {
            chunk->set_material(l_pos.x, l_pos.y, l_pos.z, (uint8_t)material_index);
        }
    }
}

int DensityGrid::get_material_id(const Vector3i &pos) const {
    int index = get_index(pos);
    if (index != -1 && index < world_material_grid.size()) {
        return (int)world_material_grid[index];
    }
    return 0;
}

void DensityGrid::set_world_material_grid(const PackedByteArray &p_grid) {
    world_material_grid = p_grid;
}

PackedByteArray DensityGrid::get_world_material_grid() const {
    return world_material_grid;
}

void DensityGrid::set_grid_size_x(int p_x) { grid_size_x = p_x > 0 ? p_x : 0; }
int DensityGrid::get_grid_size_x() const { return grid_size_x; }

void DensityGrid::set_grid_size_y(int p_y) { grid_size_y = p_y > 0 ? p_y : 0; }
int DensityGrid::get_grid_size_y() const { return grid_size_y; }

void DensityGrid::set_grid_size_z(int p_z) { grid_size_z = p_z > 0 ? p_z : 0; }
int DensityGrid::get_grid_size_z() const { return grid_size_z; }

void DensityGrid::set_surface_threshold(float p_threshold) { surface_threshold = p_threshold; }
float DensityGrid::get_surface_threshold() const { return surface_threshold; }

void DensityGrid::set_zone_at(const Vector3i &pos, int zone_id) {
    int index = get_index(pos);
    if (index != -1 && index < world_zone_grid.size()) {
        world_zone_grid[index] = zone_id;

        Vector3i l_pos = world_to_local_coord(pos);
        VoxelChunk *chunk = get_chunk_for_voxel(pos, true);
        if (chunk) {
            chunk->set_zone(l_pos.x, l_pos.y, l_pos.z, zone_id);
        }
    }
}

int DensityGrid::get_zone_at(const Vector3i &pos) const {
    int index = get_index(pos);
    if (index != -1 && index < world_zone_grid.size()) {
        return world_zone_grid[index];
    }
    return 0;
}

String DensityGrid::get_zone_name_by_id(int zone_id) const {
    if (zone_id >= 0 && zone_id < (int)zone_id_map.size()) {
        return zone_id_map[zone_id];
    }
    return "undefined";
}

int DensityGrid::get_zone_count() const {
    return (int)zone_id_map.size();
}

int DensityGrid::register_zone_name(String name) {
    for (size_t i = 0; i < zone_id_map.size(); ++i) {
        if (zone_id_map[i] == name) {
            return (int)i;
        }
    }
    zone_id_map.push_back(name);
    return (int)zone_id_map.size() - 1;
}

void DensityGrid::set_hermite_edge(const Vector3i &pos, int axis, float t, const Vector3 &normal) {
    uint64_t key = get_hermite_edge_key(pos.x, pos.y, pos.z, axis);
    HermiteEdgeData data;
    data.t = Math::clamp(t, 0.0f, 1.0f);
    data.normal = normal.length_squared() > 1e-8f ? normal.normalized() : Vector3(0, 1, 0);
    hermite_edges[key] = data;
}

bool DensityGrid::get_hermite_edge(const Vector3i &pos, int axis, HermiteEdgeData &out_data) const {
    uint64_t key = get_hermite_edge_key(pos.x, pos.y, pos.z, axis);
    auto it = hermite_edges.find(key);
    if (it != hermite_edges.end()) {
        out_data = it->second;
        return true;
    }
    return false;
}

void DensityGrid::clear_hermite_edges() {
    hermite_edges.clear();
}

} // namespace godot
