#include "undergen_scene_spawner_node.h"
#include <godot_cpp/variant/utility_functions.hpp>
#include <godot_cpp/core/math.hpp>
#include <godot_cpp/classes/resource_loader.hpp>

namespace godot {

UnderGenSceneSpawnerNode::UnderGenSceneSpawnerNode() {
    rng.instantiate();
}
UnderGenSceneSpawnerNode::~UnderGenSceneSpawnerNode() {}

void UnderGenSceneSpawnerNode::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_scene_to_spawn", "scene"), &UnderGenSceneSpawnerNode::set_scene_to_spawn);
    ClassDB::bind_method(D_METHOD("get_scene_to_spawn"), &UnderGenSceneSpawnerNode::get_scene_to_spawn);
    ClassDB::bind_method(D_METHOD("set_spawn_probability", "prob"), &UnderGenSceneSpawnerNode::set_spawn_probability);
    ClassDB::bind_method(D_METHOD("get_spawn_probability"), &UnderGenSceneSpawnerNode::get_spawn_probability);
    ClassDB::bind_method(D_METHOD("set_random_y_rotation", "enabled"), &UnderGenSceneSpawnerNode::set_random_y_rotation);
    ClassDB::bind_method(D_METHOD("get_random_y_rotation"), &UnderGenSceneSpawnerNode::get_random_y_rotation);
    ClassDB::bind_method(D_METHOD("set_random_seed", "seed"), &UnderGenSceneSpawnerNode::set_random_seed);
    ClassDB::bind_method(D_METHOD("get_random_seed"), &UnderGenSceneSpawnerNode::get_random_seed);

    ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "scene_to_spawn", PROPERTY_HINT_RESOURCE_TYPE, "PackedScene"), "set_scene_to_spawn", "get_scene_to_spawn");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "spawn_probability", PROPERTY_HINT_RANGE, "0,1,0.01"), "set_spawn_probability", "get_spawn_probability");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "random_y_rotation"), "set_random_y_rotation", "get_random_y_rotation");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "random_seed"), "set_random_seed", "get_random_seed");
}

void UnderGenSceneSpawnerNode::set_scene_to_spawn(const Variant &p_scene) {
    if (p_scene.get_type() == Variant::STRING) {
        String path = p_scene;
        if (!path.is_empty()) {
            scene_to_spawn = ResourceLoader::get_singleton()->load(path);
        } else {
            scene_to_spawn = Ref<PackedScene>();
        }
    } else if (p_scene.get_type() == Variant::OBJECT) {
        scene_to_spawn = p_scene;
    } else {
        scene_to_spawn = Ref<PackedScene>();
    }
}
Ref<PackedScene> UnderGenSceneSpawnerNode::get_scene_to_spawn() const { return scene_to_spawn; }
void UnderGenSceneSpawnerNode::set_spawn_probability(float p_prob) { spawn_probability = Math::clamp(p_prob, 0.0f, 1.0f); }
float UnderGenSceneSpawnerNode::get_spawn_probability() const { return spawn_probability; }
void UnderGenSceneSpawnerNode::set_random_y_rotation(bool p_enabled) { random_y_rotation = p_enabled; }
bool UnderGenSceneSpawnerNode::get_random_y_rotation() const { return random_y_rotation; }
void UnderGenSceneSpawnerNode::set_random_seed(int64_t p_seed) { random_seed = p_seed; rng->set_seed(p_seed); }
int64_t UnderGenSceneSpawnerNode::get_random_seed() const { return random_seed; }

void UnderGenSceneSpawnerNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    // Passthrough - actual spawning requires scene tree, see execute_with_parent
    outputs[0] = inputs.get(0, Ref<UnderGenPointSet>());
}

void UnderGenSceneSpawnerNode::execute_with_parent(const Dictionary &inputs, Dictionary &outputs, Node3D* parent_node) {
    if (!parent_node) {
        UtilityFunctions::printerr("UnderGenSceneSpawnerNode: No parent node.");
        return;
    }

    Ref<UnderGenPointSet> point_set = inputs.get(0, Ref<UnderGenPointSet>());
    if (point_set.is_null()) {
        UtilityFunctions::printerr("UnderGenSceneSpawnerNode: No PointSet on port 0.");
        return;
    }

    rng->set_seed(random_seed);
    int spawned_count = 0;

    const auto& raw = point_set->get_points_raw();
    for (const auto& p : raw) {
        // Probabilistic rejection
        if (rng->randf() > spawn_probability * p.density) continue;

        // Single scene: every point gets the same PackedScene
        if (scene_to_spawn.is_null()) {
            UtilityFunctions::printerr("UnderGenSceneSpawnerNode: No scene_to_spawn assigned.");
            break;
        }

        Node* instance = scene_to_spawn->instantiate();
        if (!instance) continue;

        Node3D* node3d = Object::cast_to<Node3D>(instance);
        if (!node3d) {
            instance->queue_free();
            continue;
        }

        Transform3D xform = p.transform;
        if (random_y_rotation) {
            float angle = rng->randf() * Math_TAU;
            xform.basis = xform.basis.rotated(Vector3(0, 1, 0), angle);
        }
        node3d->set_transform(xform);

        parent_node->add_child(node3d);
        spawned_count++;
    }

    UtilityFunctions::print("UnderGenSceneSpawnerNode: Spawned ", spawned_count, " entities.");
    outputs[0] = point_set; // Pass through for chaining
}

} // namespace godot
