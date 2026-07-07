#ifndef UNDERGEN_MESHER_NODE_H
#define UNDERGEN_MESHER_NODE_H

#include "undergen_node.h"
#include <godot_cpp/variant/typed_array.hpp>
#include <godot_cpp/classes/material.hpp>
#include <godot_cpp/classes/rd_shader_file.hpp>
#include <godot_cpp/classes/node3d.hpp>

namespace godot {

// Terminal node: reads the Generation Context and spawns MCChunk nodes
// as children of a provided parent Node3D.
class UnderGenMesherNode : public UnderGenNode {
    GDCLASS(UnderGenMesherNode, UnderGenNode);

private:
    int chunk_size = 16;
    float voxel_size = 1.0f;
    bool generate_collision = true;
    bool generate_occluder = false;
    TypedArray<Material> terrain_materials;
    Ref<Material> liquid_material;
    int liquid_material_id = 2;
    bool generate_liquid_trigger = true;
    Ref<RDShaderFile> compute_shader;

protected:
    static void _bind_methods();

public:
    UnderGenMesherNode();
    virtual ~UnderGenMesherNode();

    void set_chunk_size(int p_size);
    int get_chunk_size() const;
    void set_voxel_size(float p_size);
    float get_voxel_size() const;
    void set_generate_collision(bool p_enabled);
    bool get_generate_collision() const;
    void set_generate_occluder(bool p_enabled);
    bool get_generate_occluder() const;
    void set_terrain_materials(const TypedArray<Material> &p_materials);
    TypedArray<Material> get_terrain_materials() const;
    void set_liquid_material(const Ref<Material> &p_material);
    Ref<Material> get_liquid_material() const;
    void set_liquid_material_id(int p_id);
    int get_liquid_material_id() const;
    void set_generate_liquid_trigger(bool p_enabled);
    bool get_generate_liquid_trigger() const;
    void set_compute_shader(const Ref<RDShaderFile> &p_shader);
    Ref<RDShaderFile> get_compute_shader() const;

    // Called with the world node as context owner
    void execute_with_parent(const Dictionary &inputs, Dictionary &outputs, Node3D* parent_node);
    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
};

} // namespace godot

#endif // UNDERGEN_MESHER_NODE_H
