// src/dc_chunk.h
#ifndef DC_CHUNK_H
#define DC_CHUNK_H

#include <godot_cpp/classes/mesh_instance3d.hpp>
#include <godot_cpp/classes/array_mesh.hpp>
#include <godot_cpp/classes/material.hpp>
#include <godot_cpp/variant/vector3i.hpp>
#include <godot_cpp/variant/packed_vector3_array.hpp>
#include <godot_cpp/variant/packed_int32_array.hpp>
#include <godot_cpp/variant/packed_vector2_array.hpp>
#include <godot_cpp/variant/dictionary.hpp>
#include <godot_cpp/core/object.hpp>
#include <godot_cpp/classes/concave_polygon_shape3d.hpp>
#include <godot_cpp/classes/collision_shape3d.hpp>
#include <godot_cpp/classes/occluder_instance3d.hpp>
#include <godot_cpp/classes/array_occluder3d.hpp>
#include <godot_cpp/classes/rd_shader_file.hpp>
#include <godot_cpp/classes/rendering_device.hpp>
#include <godot_cpp/classes/rendering_server.hpp>

#include "density_grid.h"

namespace godot {

class DCChunk : public MeshInstance3D {
    GDCLASS(DCChunk, MeshInstance3D);

private:
    int chunk_size = 16;
    float voxel_size = 1.0;
    Vector3i chunk_grid_offset = Vector3i(0, 0, 0);
    TypedArray<Material> materials;
    bool smooth_normals = false;
    bool flip_normals = false;
    bool generate_collision = true;
    bool generate_occluder = true;
    Ref<RDShaderFile> compute_shader; // Added for interface compatibility

    Ref<DensityGrid> density_grid;
    Ref<ConcavePolygonShape3D> collision_shape;
    Ref<ArrayOccluder3D> occluder_shape;

    // Liquid Properties
    int liquid_material_id = 2;
    int flow_spread_limit = 7;
    Ref<Material> liquid_material;
    bool generate_liquid_trigger = true;

    // Dual Contouring Specific
    bool use_qef = true;
    float qef_regularization = 1e-4f;

    void _generate_liquid_mesh();
    void _clear_liquid();

    // --- Private method declarations for collision ---
    void _generate_collision(const PackedVector3Array &p_vertices, const PackedInt32Array &p_indices);
    void _clear_collision();

    // --- Private method declarations for occluder ---
    void _generate_occluder(const PackedVector3Array &p_vertices, const PackedInt32Array &p_indices);
    void _clear_occluder();

protected:
    static void _bind_methods();
    void _notification(int p_what);

public:
    DCChunk();
    ~DCChunk();

    void generate_mesh_from_density_grid();

    void set_chunk_size(int p_size);
    int get_chunk_size() const;

    void set_voxel_size(float p_size);
    float get_voxel_size() const;

    void set_chunk_grid_offset(const Vector3i &p_offset);
    Vector3i get_chunk_grid_offset() const;

    void set_density_grid(const Ref<DensityGrid> &p_grid);
    Ref<DensityGrid> get_density_grid() const;

    void set_generate_collision(bool p_generate);
    bool get_generate_collision() const;

    void set_generate_occluder(bool p_generate);
    bool get_generate_occluder() const;

    void set_materials(const TypedArray<Material> &p_materials);
    TypedArray<Material> get_materials() const;

    void set_compute_shader(const Ref<RDShaderFile> &p_shader);
    Ref<RDShaderFile> get_compute_shader() const;

    void set_liquid_material(const Ref<Material> &p_material);
    Ref<Material> get_liquid_material() const;
    void set_liquid_material_id(int p_id);
    int get_liquid_material_id() const;
    void set_generate_liquid_trigger(bool p_enabled);
    bool get_generate_liquid_trigger() const;
    void set_flow_spread_limit(int p_limit);
    int get_flow_spread_limit() const;

    void set_smooth_normals(bool p_smooth);
    bool get_smooth_normals() const;

    void set_flip_normals(bool p_flip);
    bool get_flip_normals() const;

    // Dual Contouring parameters
    void set_use_qef(bool p_use);
    bool get_use_qef() const;

    void set_qef_regularization(float p_reg);
    float get_qef_regularization() const;
};

} // namespace godot

#endif // DC_CHUNK_H
