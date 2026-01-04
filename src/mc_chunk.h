// src/mc_chunk.h
#ifndef MC_CHUNK_H
#define MC_CHUNK_H

#include <godot_cpp/classes/mesh_instance3d.hpp>
#include <godot_cpp/classes/array_mesh.hpp>
#include <godot_cpp/classes/material.hpp>
#include <godot_cpp/variant/vector3i.hpp>
#include <godot_cpp/variant/packed_vector3_array.hpp>
#include <godot_cpp/variant/packed_int32_array.hpp>
#include <godot_cpp/variant/packed_vector2_array.hpp>
#include <godot_cpp/variant/dictionary.hpp>
#include <godot_cpp/core/object.hpp>
#include <godot_cpp/classes/concave_polygon_shape3d.hpp> // Added for collision shape
#include <godot_cpp/classes/collision_shape3d.hpp>     // Added for collision node
#include <godot_cpp/classes/occluder_instance3d.hpp>    // For occluder node
#include <godot_cpp/classes/array_occluder3d.hpp>     // For occluder shape
#include <godot_cpp/classes/rd_shader_file.hpp>
#include <godot_cpp/classes/rendering_device.hpp>
#include <godot_cpp/classes/rendering_server.hpp>

#include "density_grid.h"
#include "level_density_grid.h"
#include "mc_tables.h"

namespace godot {

class MCChunk : public MeshInstance3D {
    GDCLASS(MCChunk, MeshInstance3D);

private:
    int chunk_size = 16;
    float voxel_size = 1.0;
    Vector3i chunk_grid_offset = Vector3i(0, 0, 0);
    Ref<Material> terrain_material;
    bool smooth_normals = true;
    bool generate_collision = true;
    bool generate_occluder = true; // Added property for occluder generation
    Ref<RDShaderFile> compute_shader;

    Ref<DensityGrid> density_grid;
    Ref<ConcavePolygonShape3D> collision_shape; // Member to hold the collision shape resource
    Ref<ArrayOccluder3D> occluder_shape;      // Member to hold the occluder shape resource

    Dictionary _march_cubes();
    Vector3 _interpolate_vertex(const Vector3 &p1, const Vector3 &p2, float val1, float val2);
    void _generate_mesh_with_compute();

    // --- Private method declarations for collision ---
    void _generate_collision(const PackedVector3Array &p_vertices, const PackedInt32Array &p_indices);
    void _clear_collision(); // Method to remove existing collision nodes

    // --- Private method declarations for occluder ---
    void _generate_occluder(const PackedVector3Array &p_vertices, const PackedInt32Array &p_indices);
    void _clear_occluder(); // Method to remove existing occluder nodes

protected:
    static void _bind_methods();
    void _notification(int p_what);

public:
    MCChunk();
    ~MCChunk();

    void generate_mesh_from_density_grid();

    void set_chunk_size(int p_size);
    int get_chunk_size() const;

    void set_voxel_size(float p_size);
    float get_voxel_size() const;

    void set_chunk_grid_offset(const Vector3i &p_offset);
    Vector3i get_chunk_grid_offset() const;

    void set_terrain_material(const Ref<Material> &p_material);
    Ref<Material> get_terrain_material() const;

    void set_density_grid(const Ref<DensityGrid> &p_grid);
    Ref<DensityGrid> get_density_grid() const;

    // New getter/setter for generate_collision
    void set_generate_collision(bool p_generate);
    bool get_generate_collision() const;

    // New getter/setter for generate_occluder
    void set_generate_occluder(bool p_generate);
    bool get_generate_occluder() const;

    // Compute Shader
    void set_compute_shader(const Ref<RDShaderFile> &p_shader);
    Ref<RDShaderFile> get_compute_shader() const;
};

} // namespace godot

#endif // MC_CHUNK_H
