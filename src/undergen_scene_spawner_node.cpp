#include "undergen_scene_spawner_node.h"
#include <godot_cpp/variant/utility_functions.hpp>
#include <godot_cpp/core/math.hpp>
#include <godot_cpp/classes/resource_loader.hpp>
#include <algorithm>

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
    ClassDB::bind_method(D_METHOD("set_align_with_normal", "align"), &UnderGenSceneSpawnerNode::set_align_with_normal);
    ClassDB::bind_method(D_METHOD("get_align_with_normal"), &UnderGenSceneSpawnerNode::get_align_with_normal);
    ClassDB::bind_method(D_METHOD("set_spawn_limit", "limit"), &UnderGenSceneSpawnerNode::set_spawn_limit);
    ClassDB::bind_method(D_METHOD("get_spawn_limit"), &UnderGenSceneSpawnerNode::get_spawn_limit);
    ClassDB::bind_method(D_METHOD("set_shuffle_points", "shuffle"), &UnderGenSceneSpawnerNode::set_shuffle_points);
    ClassDB::bind_method(D_METHOD("get_shuffle_points"), &UnderGenSceneSpawnerNode::get_shuffle_points);
    ClassDB::bind_method(D_METHOD("set_random_seed", "seed"), &UnderGenSceneSpawnerNode::set_random_seed);
    ClassDB::bind_method(D_METHOD("get_random_seed"), &UnderGenSceneSpawnerNode::get_random_seed);
    ClassDB::bind_method(D_METHOD("set_consume_points", "consume"), &UnderGenSceneSpawnerNode::set_consume_points);
    ClassDB::bind_method(D_METHOD("get_consume_points"), &UnderGenSceneSpawnerNode::get_consume_points);
    ClassDB::bind_method(D_METHOD("set_parent_node", "parent_node"), &UnderGenSceneSpawnerNode::set_parent_node);
    ClassDB::bind_method(D_METHOD("get_parent_node"), &UnderGenSceneSpawnerNode::get_parent_node);
    ClassDB::bind_method(D_METHOD("set_force_readable_names", "force"), &UnderGenSceneSpawnerNode::set_force_readable_names);
    ClassDB::bind_method(D_METHOD("get_force_readable_names"), &UnderGenSceneSpawnerNode::get_force_readable_names);
    ClassDB::bind_method(D_METHOD("set_multiplayer_spawner", "path"), &UnderGenSceneSpawnerNode::set_multiplayer_spawner);
    ClassDB::bind_method(D_METHOD("get_multiplayer_spawner"), &UnderGenSceneSpawnerNode::get_multiplayer_spawner);
    ClassDB::bind_method(D_METHOD("spawn_scenes_with_parent", "inputs", "parent_node"), &UnderGenSceneSpawnerNode::spawn_scenes_with_parent);

    ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "scene_to_spawn", PROPERTY_HINT_RESOURCE_TYPE, "PackedScene"), "set_scene_to_spawn", "get_scene_to_spawn");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "spawn_probability", PROPERTY_HINT_RANGE, "0,1,0.01"), "set_spawn_probability", "get_spawn_probability");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "random_y_rotation"), "set_random_y_rotation", "get_random_y_rotation");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "align_with_normal"), "set_align_with_normal", "get_align_with_normal");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "spawn_limit"), "set_spawn_limit", "get_spawn_limit");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "shuffle_points"), "set_shuffle_points", "get_shuffle_points");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "random_seed"), "set_random_seed", "get_random_seed");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "consume_points"), "set_consume_points", "get_consume_points");
    ADD_PROPERTY(PropertyInfo(Variant::NODE_PATH, "parent_node"), "set_parent_node", "get_parent_node");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "force_readable_names"), "set_force_readable_names", "get_force_readable_names");
    ADD_PROPERTY(PropertyInfo(Variant::NODE_PATH, "multiplayer_spawner"), "set_multiplayer_spawner", "get_multiplayer_spawner");
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
void UnderGenSceneSpawnerNode::set_align_with_normal(bool p_align) { align_with_normal = p_align; }
bool UnderGenSceneSpawnerNode::get_align_with_normal() const { return align_with_normal; }
void UnderGenSceneSpawnerNode::set_spawn_limit(int p_limit) { spawn_limit = p_limit < 0 ? 0 : p_limit; }
int UnderGenSceneSpawnerNode::get_spawn_limit() const { return spawn_limit; }
void UnderGenSceneSpawnerNode::set_shuffle_points(bool p_shuffle) { shuffle_points = p_shuffle; }
bool UnderGenSceneSpawnerNode::get_shuffle_points() const { return shuffle_points; }
void UnderGenSceneSpawnerNode::set_random_seed(int64_t p_seed) { random_seed = p_seed; rng->set_seed(p_seed); }
int64_t UnderGenSceneSpawnerNode::get_random_seed() const { return random_seed; }
void UnderGenSceneSpawnerNode::set_consume_points(bool p_consume) { consume_points = p_consume; }
bool UnderGenSceneSpawnerNode::get_consume_points() const { return consume_points; }
void UnderGenSceneSpawnerNode::set_parent_node(const NodePath &p_path) { parent_node_path = p_path; }
NodePath UnderGenSceneSpawnerNode::get_parent_node() const { return parent_node_path; }
void UnderGenSceneSpawnerNode::set_force_readable_names(bool p_force) { force_readable_names = p_force; }
bool UnderGenSceneSpawnerNode::get_force_readable_names() const { return force_readable_names; }
void UnderGenSceneSpawnerNode::set_multiplayer_spawner(const NodePath &p_path) { multiplayer_spawner_path = p_path; }
NodePath UnderGenSceneSpawnerNode::get_multiplayer_spawner() const { return multiplayer_spawner_path; }

void UnderGenSceneSpawnerNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    // Passthrough - actual spawning requires scene tree, see execute_with_parent
    outputs[0] = inputs.get(0, Ref<UnderGenPointSet>());
}

void UnderGenSceneSpawnerNode::execute_with_parent(const Dictionary &inputs, Dictionary &outputs, Node3D* parent_node) {
    if (!parent_node) {
        UtilityFunctions::printerr("UnderGenSceneSpawnerNode: No parent node.");
        return;
    }

    Node *target_parent = parent_node;
    if (!parent_node_path.is_empty()) {
        Node *custom_parent = parent_node->get_node_or_null(parent_node_path);
        if (custom_parent) {
            target_parent = custom_parent;
        } else {
            UtilityFunctions::printerr("UnderGenSceneSpawnerNode: Custom parent node not found at path: ", parent_node_path);
        }
    }

    MultiplayerSpawner *spawner = nullptr;
    if (!multiplayer_spawner_path.is_empty()) {
        Node *spawner_node = parent_node->get_node_or_null(multiplayer_spawner_path);
        if (spawner_node) {
            spawner = Object::cast_to<MultiplayerSpawner>(spawner_node);
            if (!spawner) {
                UtilityFunctions::printerr("UnderGenSceneSpawnerNode: Node at multiplayer_spawner path is not a MultiplayerSpawner.");
            }
        } else {
            UtilityFunctions::printerr("UnderGenSceneSpawnerNode: MultiplayerSpawner not found at path: ", multiplayer_spawner_path);
        }
    }

    Ref<UnderGenPointSet> point_set = inputs.get(0, Ref<UnderGenPointSet>());
    if (point_set.is_null()) {
        UtilityFunctions::printerr("UnderGenSceneSpawnerNode: No PointSet on port 0.");
        return;
    }

    rng->set_seed(random_seed);
    int spawned_count = 0;

    const auto& raw = point_set->get_points_raw();
    if (raw.empty()) {
        outputs[0] = point_set;
        return;
    }

    std::vector<size_t> indices(raw.size());
    for (size_t i = 0; i < raw.size(); ++i) {
        indices[i] = i;
    }

    if (shuffle_points) {
        for (size_t i = raw.size() - 1; i > 0; --i) {
            size_t j = rng->randi() % (i + 1);
            std::swap(indices[i], indices[j]);
        }
    }

    std::vector<size_t> spawned_indices;

    for (size_t idx : indices) {
        if (spawn_limit > 0 && spawned_count >= spawn_limit) {
            break;
        }

        const auto& p = raw[idx];

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
        if (align_with_normal) {
            Vector3 normal = p.attributes.get("normal", Vector3(0, 1, 0));
            Vector3 up = Vector3(0, 1, 0);
            if (normal.is_equal_approx(-up)) {
                xform.basis = xform.basis.rotated(Vector3(1, 0, 0), Math_PI);
            } else if (!normal.is_equal_approx(up)) {
                Vector3 axis = up.cross(normal).normalized();
                float angle = Math::acos(up.dot(normal));
                xform.basis = xform.basis.rotated(axis, angle);
            }

            if (random_y_rotation) {
                float angle = rng->randf() * Math_TAU;
                xform.basis = xform.basis.rotated(normal, angle);
            }
        } else {
            if (random_y_rotation) {
                float angle = rng->randf() * Math_TAU;
                xform.basis = xform.basis.rotated(Vector3(0, 1, 0), angle);
            }
        }
        if (spawner) {
            Dictionary spawn_data;
            spawn_data["scene_path"] = scene_to_spawn->get_path();
            spawn_data["transform"] = xform;

            Node *spawned = spawner->spawn(spawn_data);
            if (spawned) {
                spawned_count++;
                spawned_indices.push_back(idx);
            } else {
                UtilityFunctions::printerr("UnderGenSceneSpawnerNode: MultiplayerSpawner failed to spawn node.");
            }
            node3d->queue_free();
        } else {
            node3d->set_transform(xform);
            target_parent->add_child(node3d, force_readable_names);
            spawned_count++;
            spawned_indices.push_back(idx);
        }
    }

    if (consume_points && !spawned_indices.empty()) {
        std::sort(spawned_indices.begin(), spawned_indices.end(), std::greater<size_t>());
        for (size_t idx : spawned_indices) {
            point_set->remove_point((int)idx);
        }
    }

    UtilityFunctions::print("UnderGenSceneSpawnerNode: Spawned ", spawned_count, " entities.");
    outputs[0] = point_set; // Pass through for chaining
}

void UnderGenSceneSpawnerNode::spawn_scenes_with_parent(const Dictionary &inputs, Node *parent_node) {
    if (!parent_node) {
        UtilityFunctions::printerr("UnderGenSceneSpawnerNode::spawn_scenes_with_parent: parent_node is null.");
        return;
    }
    Node3D *parent_3d = Object::cast_to<Node3D>(parent_node);
    if (!parent_3d) {
        UtilityFunctions::printerr("UnderGenSceneSpawnerNode::spawn_scenes_with_parent: parent_node is not a Node3D.");
        return;
    }
    Dictionary dummy_outputs;
    execute_with_parent(inputs, dummy_outputs, parent_3d);
}

} // namespace godot
