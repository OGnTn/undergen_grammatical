#ifndef UNDERGEN_SCENE_SPAWNER_NODE_H
#define UNDERGEN_SCENE_SPAWNER_NODE_H

#include "undergen_node.h"
#include "undergen_point_set.h"
#include <godot_cpp/classes/packed_scene.hpp>
#include <godot_cpp/classes/random_number_generator.hpp>
#include <godot_cpp/classes/node3d.hpp>

namespace godot {

// Terminal node: spawns PackedScene instances at positions given by a PointSet.
// Respects density as a spawn probability (0=never, 1=always).
// Supports zone→scene mapping: if a point's zone_name matches a key in zone_scene_map,
// that scene is used instead of scene_to_spawn.
class UnderGenSceneSpawnerNode : public UnderGenNode {
    GDCLASS(UnderGenSceneSpawnerNode, UnderGenNode);

private:
    Ref<PackedScene> scene_to_spawn; // Default / fallback scene
    Dictionary zone_scene_map;       // { "zone_name": PackedScene, ... }
    float spawn_probability = 1.0f;
    bool random_y_rotation = true;
    int64_t random_seed = 0;
    Ref<RandomNumberGenerator> rng;

protected:
    static void _bind_methods();

public:
    UnderGenSceneSpawnerNode();
    virtual ~UnderGenSceneSpawnerNode();

    void set_scene_to_spawn(const Ref<PackedScene> &p_scene);
    Ref<PackedScene> get_scene_to_spawn() const;
    void set_zone_scene_map(const Dictionary &p_map);
    Dictionary get_zone_scene_map() const;
    void set_spawn_probability(float p_prob);
    float get_spawn_probability() const;
    void set_random_y_rotation(bool p_enabled);
    bool get_random_y_rotation() const;
    void set_random_seed(int64_t p_seed);
    int64_t get_random_seed() const;

    void execute_with_parent(const Dictionary &inputs, Dictionary &outputs, Node3D* parent_node);
    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
};

} // namespace godot

#endif // UNDERGEN_SCENE_SPAWNER_NODE_H
