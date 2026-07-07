#include "undergen_mesher_node.h"
#include "mc_chunk.h"
#include "density_grid.h"
#include <godot_cpp/variant/utility_functions.hpp>
#include <godot_cpp/core/math.hpp>

namespace godot {

UnderGenMesherNode::UnderGenMesherNode() {}
UnderGenMesherNode::~UnderGenMesherNode() {}

void UnderGenMesherNode::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_chunk_size", "size"), &UnderGenMesherNode::set_chunk_size);
    ClassDB::bind_method(D_METHOD("get_chunk_size"), &UnderGenMesherNode::get_chunk_size);
    ClassDB::bind_method(D_METHOD("set_voxel_size", "size"), &UnderGenMesherNode::set_voxel_size);
    ClassDB::bind_method(D_METHOD("get_voxel_size"), &UnderGenMesherNode::get_voxel_size);
    ClassDB::bind_method(D_METHOD("set_generate_collision", "enabled"), &UnderGenMesherNode::set_generate_collision);
    ClassDB::bind_method(D_METHOD("get_generate_collision"), &UnderGenMesherNode::get_generate_collision);
    ClassDB::bind_method(D_METHOD("set_generate_occluder", "enabled"), &UnderGenMesherNode::set_generate_occluder);
    ClassDB::bind_method(D_METHOD("get_generate_occluder"), &UnderGenMesherNode::get_generate_occluder);
    ClassDB::bind_method(D_METHOD("set_terrain_materials", "materials"), &UnderGenMesherNode::set_terrain_materials);
    ClassDB::bind_method(D_METHOD("get_terrain_materials"), &UnderGenMesherNode::get_terrain_materials);
    ClassDB::bind_method(D_METHOD("set_liquid_material", "material"), &UnderGenMesherNode::set_liquid_material);
    ClassDB::bind_method(D_METHOD("get_liquid_material"), &UnderGenMesherNode::get_liquid_material);
    ClassDB::bind_method(D_METHOD("set_liquid_material_id", "id"), &UnderGenMesherNode::set_liquid_material_id);
    ClassDB::bind_method(D_METHOD("get_liquid_material_id"), &UnderGenMesherNode::get_liquid_material_id);
    ClassDB::bind_method(D_METHOD("set_generate_liquid_trigger", "enabled"), &UnderGenMesherNode::set_generate_liquid_trigger);
    ClassDB::bind_method(D_METHOD("get_generate_liquid_trigger"), &UnderGenMesherNode::get_generate_liquid_trigger);
    ClassDB::bind_method(D_METHOD("set_compute_shader", "shader"), &UnderGenMesherNode::set_compute_shader);
    ClassDB::bind_method(D_METHOD("get_compute_shader"), &UnderGenMesherNode::get_compute_shader);

    ADD_PROPERTY(PropertyInfo(Variant::INT, "chunk_size"), "set_chunk_size", "get_chunk_size");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "voxel_size"), "set_voxel_size", "get_voxel_size");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "generate_collision"), "set_generate_collision", "get_generate_collision");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "generate_occluder"), "set_generate_occluder", "get_generate_occluder");
    ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "terrain_materials"), "set_terrain_materials", "get_terrain_materials");
    ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "liquid_material", PROPERTY_HINT_RESOURCE_TYPE, "Material"), "set_liquid_material", "get_liquid_material");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "liquid_material_id", PROPERTY_HINT_RANGE, "0,255,1"), "set_liquid_material_id", "get_liquid_material_id");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "generate_liquid_trigger"), "set_generate_liquid_trigger", "get_generate_liquid_trigger");
    ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "compute_shader", PROPERTY_HINT_RESOURCE_TYPE, "RDShaderFile"), "set_compute_shader", "get_compute_shader");
}

void UnderGenMesherNode::set_chunk_size(int p_size) { chunk_size = p_size; }
int UnderGenMesherNode::get_chunk_size() const { return chunk_size; }
void UnderGenMesherNode::set_voxel_size(float p_size) { voxel_size = p_size; }
float UnderGenMesherNode::get_voxel_size() const { return voxel_size; }
void UnderGenMesherNode::set_generate_collision(bool p_enabled) { generate_collision = p_enabled; }
bool UnderGenMesherNode::get_generate_collision() const { return generate_collision; }
void UnderGenMesherNode::set_generate_occluder(bool p_enabled) { generate_occluder = p_enabled; }
bool UnderGenMesherNode::get_generate_occluder() const { return generate_occluder; }
void UnderGenMesherNode::set_terrain_materials(const TypedArray<Material> &p_materials) { terrain_materials = p_materials; }
TypedArray<Material> UnderGenMesherNode::get_terrain_materials() const { return terrain_materials; }
void UnderGenMesherNode::set_liquid_material(const Ref<Material> &p_material) { liquid_material = p_material; }
Ref<Material> UnderGenMesherNode::get_liquid_material() const { return liquid_material; }
void UnderGenMesherNode::set_liquid_material_id(int p_id) { liquid_material_id = Math::clamp(p_id, 0, 255); }
int UnderGenMesherNode::get_liquid_material_id() const { return liquid_material_id; }
void UnderGenMesherNode::set_generate_liquid_trigger(bool p_enabled) { generate_liquid_trigger = p_enabled; }
bool UnderGenMesherNode::get_generate_liquid_trigger() const { return generate_liquid_trigger; }
void UnderGenMesherNode::set_compute_shader(const Ref<RDShaderFile> &p_shader) { compute_shader = p_shader; }
Ref<RDShaderFile> UnderGenMesherNode::get_compute_shader() const { return compute_shader; }

void UnderGenMesherNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    // This base version stores the context. The parent UnderGenWorld3D
    // calls execute_with_parent once the scene tree is available.
    outputs[0] = inputs.get(0, Dictionary());
}

void UnderGenMesherNode::execute_with_parent(const Dictionary &inputs, Dictionary &outputs, Node3D* parent_node) {
    if (!parent_node) {
        UtilityFunctions::printerr("UnderGenMesherNode: No parent node provided.");
        return;
    }

    Dictionary context = inputs.get(0, Dictionary());
    if (context.is_empty()) return;

    Ref<DensityGrid> grid = context.get("grid", Ref<DensityGrid>());
    if (grid.is_null()) {
        UtilityFunctions::printerr("UnderGenMesherNode: Grid not found in context.");
        return;
    }

    Vector3i dims = grid->get_grid_dimensions();
    int count_x = Math::max(1, (int)Math::ceil((float)(dims.x - 1) / chunk_size));
    int count_y = Math::max(1, (int)Math::ceil((float)(dims.y - 1) / chunk_size));
    int count_z = Math::max(1, (int)Math::ceil((float)(dims.z - 1) / chunk_size));

    UtilityFunctions::print("UnderGenMesherNode: Spawning ", count_x * count_y * count_z, " terrain chunks...");

    for (int x = 0; x < count_x; ++x) {
        for (int y = 0; y < count_y; ++y) {
            for (int z = 0; z < count_z; ++z) {
                MCChunk* chunk = memnew(MCChunk);
                chunk->set_name("TerrainChunk_" + String::num_int64(x) + "_" + String::num_int64(y) + "_" + String::num_int64(z));
                chunk->set_chunk_size(chunk_size);
                chunk->set_voxel_size(voxel_size);
                chunk->set_chunk_grid_offset(Vector3i(x, y, z) * chunk_size);
                chunk->set_density_grid(grid);
                chunk->set_generate_collision(generate_collision);
                chunk->set_generate_occluder(generate_occluder);
                if (liquid_material.is_valid()) {
                    // We'll call set_liquid_material, set_liquid_material_id, set_generate_liquid_trigger on the chunk
                    chunk->call("set_liquid_material", liquid_material);
                    chunk->call("set_liquid_material_id", liquid_material_id);
                    chunk->call("set_generate_liquid_trigger", generate_liquid_trigger);
                }
                if (!terrain_materials.is_empty()) {
                    chunk->set_materials(terrain_materials);
                }
                if (compute_shader.is_valid()) {
                    chunk->set_compute_shader(compute_shader);
                }
                chunk->set_position(Vector3(x, y, z) * chunk_size * voxel_size);
                parent_node->add_child(chunk);
                chunk->call_deferred("generate_mesh_from_density_grid");
            }
        }
    }

    outputs[0] = context;
    UtilityFunctions::print("UnderGenMesherNode: Chunk spawning complete.");
}

} // namespace godot
