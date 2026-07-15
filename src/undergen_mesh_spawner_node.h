#ifndef UNDERGEN_MESH_SPAWNER_NODE_H
#define UNDERGEN_MESH_SPAWNER_NODE_H

#include "undergen_node.h"
#include "undergen_point_set.h"
#include <godot_cpp/classes/mesh.hpp>
#include <godot_cpp/classes/material.hpp>
#include <godot_cpp/classes/random_number_generator.hpp>
#include <godot_cpp/classes/node.hpp>
#include <godot_cpp/classes/geometry_instance3d.hpp>

namespace godot {

class UnderGenMeshSpawnerNode : public UnderGenNode {
    GDCLASS(UnderGenMeshSpawnerNode, UnderGenNode);

private:
    Ref<Mesh> mesh_to_spawn;
    Ref<Material> material_override;
    float spawn_probability = 1.0f;
    int instances_per_point = 1;
    float max_distance_from_point = 0.0f;
    float min_y_offset = 0.0f;
    float max_y_offset = 0.0f;
    Vector3 min_rotation = Vector3(0, 0, 0); // In degrees
    Vector3 max_rotation = Vector3(0, 0, 0); // In degrees
    Vector3 min_scale = Vector3(1, 1, 1);
    Vector3 max_scale = Vector3(1, 1, 1);
    int64_t random_seed = 0;
    bool consume_points = false;
    bool align_with_normal = false;
    int cast_shadows = 1; // standard GeometryInstance3D::SHADOW_CASTING_SETTING_ON
    float chunk_size = 16.0f;

    Ref<RandomNumberGenerator> rng;

protected:
    static void _bind_methods();

public:
    UnderGenMeshSpawnerNode();
    virtual ~UnderGenMeshSpawnerNode();

    void set_mesh_to_spawn(const Ref<Mesh> &p_mesh);
    Ref<Mesh> get_mesh_to_spawn() const;

    void set_material_override(const Ref<Material> &p_material);
    Ref<Material> get_material_override() const;

    void set_spawn_probability(float p_prob);
    float get_spawn_probability() const;

    void set_instances_per_point(int p_instances);
    int get_instances_per_point() const;

    void set_max_distance_from_point(float p_dist);
    float get_max_distance_from_point() const;

    void set_min_y_offset(float p_offset);
    float get_min_y_offset() const;

    void set_max_y_offset(float p_offset);
    float get_max_y_offset() const;

    void set_min_rotation(const Vector3 &p_rot);
    Vector3 get_min_rotation() const;

    void set_max_rotation(const Vector3 &p_rot);
    Vector3 get_max_rotation() const;

    void set_min_scale(const Vector3 &p_scale);
    Vector3 get_min_scale() const;

    void set_max_scale(const Vector3 &p_scale);
    Vector3 get_max_scale() const;

    void set_random_seed(int64_t p_seed);
    int64_t get_random_seed() const;

    void set_consume_points(bool p_consume);
    bool get_consume_points() const;

    void set_align_with_normal(bool p_align);
    bool get_align_with_normal() const;

    void set_cast_shadows(int p_shadows);
    int get_cast_shadows() const;

    void set_chunk_size(float p_size);
    float get_chunk_size() const;

    void execute_with_parent(const Dictionary &inputs, Dictionary &outputs, Node *target_parent);
    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
};

} // namespace godot

#endif // UNDERGEN_MESH_SPAWNER_NODE_H
