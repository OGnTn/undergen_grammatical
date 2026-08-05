// src/density_grid.h
#ifndef DENSITY_GRID_H
#define DENSITY_GRID_H

#include <godot_cpp/classes/resource.hpp>
#include <godot_cpp/variant/packed_byte_array.hpp>
#include <godot_cpp/variant/packed_float32_array.hpp>
#include <godot_cpp/variant/packed_int32_array.hpp>
#include <godot_cpp/variant/vector3i.hpp>
#include <vector>
#include <string>
#include <unordered_map>
#include "voxel_chunk.h"

namespace godot {

struct Vector3iHash {
    std::size_t operator()(const Vector3i &v) const {
        std::size_t h1 = std::hash<int>()(v.x);
        std::size_t h2 = std::hash<int>()(v.y);
        std::size_t h3 = std::hash<int>()(v.z);
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

struct HermiteEdgeData {
    float t = 0.5f;
    Vector3 normal = Vector3(0, 1, 0);
};

inline uint64_t get_hermite_edge_key(int x, int y, int z, int axis) {
    return (((uint64_t)(x + 32768) & 0xFFFFF) << 42) |
           (((uint64_t)(y + 32768) & 0xFFFFF) << 22) |
           (((uint64_t)(z + 32768) & 0xFFFFF) << 2) |
           ((uint64_t)axis & 0x3);
}

class DensityGrid : public Resource {
    GDCLASS(DensityGrid, Resource);

protected:
    // Dynamic / Flattened cached PackedArrays for backward compatibility
    mutable PackedFloat32Array world_density_grid;
    mutable PackedByteArray world_material_grid;
    mutable PackedInt32Array world_zone_grid;
    mutable bool flat_cache_dirty = true;

    // Sparse Voxel Chunk Storage
    std::unordered_map<Vector3i, VoxelChunk*, Vector3iHash> chunks;
    std::unordered_map<uint64_t, HermiteEdgeData> hermite_edges;

    std::vector<String> zone_id_map;
    int grid_size_x = 0;
    int grid_size_y = 0;
    int grid_size_z = 0;
    float surface_threshold = 0.0f;
    float default_initial_value = 1.0f;
    int default_material_idx = 0;

protected:
    static void _bind_methods();

public:
    DensityGrid();
    ~DensityGrid();

    // Initialization & Accessors
    void initialize_grid(int size_x, int size_y, int size_z, float default_value = 1.0f, int default_material_index = 0);
    bool is_valid_position(const Vector3i &pos) const;
    int get_index(const Vector3i &pos) const;
    bool set_cell(const Vector3i &pos, float value);
    bool set_cell_fast(const Vector3i &pos, float value);
    float get_cell(const Vector3i &pos, float default_value = 1.0f) const;
    Vector3i get_grid_dimensions() const;

    // Sparse Chunk API
    static Vector3i world_to_chunk_coord(const Vector3i &pos);
    static Vector3i world_to_local_coord(const Vector3i &pos);

    VoxelChunk* get_chunk(const Vector3i &chunk_coord, bool create_if_missing = false);
    VoxelChunk* get_chunk_for_voxel(const Vector3i &pos, bool create_if_missing = false);
    const std::unordered_map<Vector3i, VoxelChunk*, Vector3iHash>& get_all_chunks() const { return chunks; }
    void clear_chunks();
    void mark_flat_cache_dirty() const { flat_cache_dirty = true; }

    // Properties exposed to Godot
    void set_world_density_grid(const PackedFloat32Array &p_grid);
    PackedFloat32Array get_world_density_grid() const;

    // C++ hot-path accessors & backward compatibility getters
    void sync_flat_cache_if_dirty() const;
    PackedFloat32Array &get_density_data_rw();
    const PackedFloat32Array &get_density_data() const;
    PackedByteArray &get_material_data_rw();
    const PackedByteArray &get_material_data() const;
    PackedInt32Array &get_zone_data_rw();
    const PackedInt32Array &get_zone_data() const;
    int64_t get_total_cell_count() const;
    int get_index_unchecked(int x, int y, int z) const;

    void set_grid_size_x(int p_x);
    int get_grid_size_x() const;

    void set_grid_size_y(int p_y);
    int get_grid_size_y() const;

    void set_grid_size_z(int p_z);
    int get_grid_size_z() const;

    void set_surface_threshold(float p_threshold);
    float get_surface_threshold() const;

    void set_material_id(const Vector3i &pos, int material_index);
    void set_material_id_fast(const Vector3i &pos, int material_index);
    int get_material_id(const Vector3i &pos) const;

    void set_world_material_grid(const PackedByteArray &p_grid);
    PackedByteArray get_world_material_grid() const;

    void set_zone_at(const Vector3i &pos, int zone_id);
    void set_zone_at_fast(const Vector3i &pos, int zone_id);
    int get_zone_at(const Vector3i &pos) const;
    String get_zone_name_by_id(int zone_id) const;
    int get_zone_count() const;
    int register_zone_name(String name);

    // Hermite Edge Data API
    void set_hermite_edge(const Vector3i &pos, int axis, float t, const Vector3 &normal);
    bool get_hermite_edge(const Vector3i &pos, int axis, HermiteEdgeData &out_data) const;
    void clear_hermite_edges();

    // SVO Raycasting API
    Dictionary raycast_svo(const Vector3 &origin, const Vector3 &dir, float max_dist, float iso_surface = 0.0f) const;
};

} // namespace godot

#endif // DENSITY_GRID_H
