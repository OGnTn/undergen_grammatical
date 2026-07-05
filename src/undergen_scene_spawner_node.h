#ifndef UNDERGEN_SCENE_SPAWNER_NODE_H
#define UNDERGEN_SCENE_SPAWNER_NODE_H

#include "undergen_node.h"
#include "undergen_point_set.h"
#include <godot_cpp/classes/packed_scene.hpp>
#include <godot_cpp/classes/random_number_generator.hpp>
#include <godot_cpp/classes/node3d.hpp>

namespace godot {

// Terminal node: spawns ONE PackedScene at every point in the input PointSet.
// For per-zone spawning, use UnderGenPointFilterNode upstream to route zones
// into separate Spawner instances, each configured with its own scene.
class UnderGenSceneSpawnerNode : public UnderGenNode {
    GDCLASS(UnderGenSceneSpawnerNode, UnderGenNode);

private:
    Ref<PackedScene> scene_to_spawn;
    float spawn_probability = 1.0f;
    bool random_y_rotation = true;
    bool align_with_normal = false;
    int spawn_limit = 0;
    bool shuffle_points = false;
    int64_t random_seed = 0;
    bool consume_points = false;
    Ref<RandomNumberGenerator> rng;

protected:
    static void _bind_methods();

public:
    UnderGenSceneSpawnerNode();
    virtual ~UnderGenSceneSpawnerNode();

    void set_scene_to_spawn(const Variant &p_scene);
    Ref<PackedScene> get_scene_to_spawn() const;
    void set_spawn_probability(float p_prob);
    float get_spawn_probability() const;
    void set_random_y_rotation(bool p_enabled);
    bool get_random_y_rotation() const;
    void set_align_with_normal(bool p_align);
    bool get_align_with_normal() const;
    void set_spawn_limit(int p_limit);
    int get_spawn_limit() const;
    void set_shuffle_points(bool p_shuffle);
    bool get_shuffle_points() const;
    void set_random_seed(int64_t p_seed);
    int64_t get_random_seed() const;
    void set_consume_points(bool p_consume);
    bool get_consume_points() const;

    void execute_with_parent(const Dictionary &inputs, Dictionary &outputs, Node3D* parent_node);
    void spawn_scenes_with_parent(const Dictionary &inputs, Node *parent_node);
    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
};

} // namespace godot

#endif // UNDERGEN_SCENE_SPAWNER_NODE_H
