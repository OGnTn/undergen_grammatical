#include "undergen_mesh_spawner_node.h"
#include <godot_cpp/variant/utility_functions.hpp>
#include <godot_cpp/core/math.hpp>
#include <godot_cpp/classes/multi_mesh.hpp>
#include <godot_cpp/classes/multi_mesh_instance3d.hpp>
#include <algorithm>
#include <tuple>
#include <map>

namespace godot {

UnderGenMeshSpawnerNode::UnderGenMeshSpawnerNode() {
    rng.instantiate();
}

UnderGenMeshSpawnerNode::~UnderGenMeshSpawnerNode() {}

void UnderGenMeshSpawnerNode::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_mesh_to_spawn", "mesh"), &UnderGenMeshSpawnerNode::set_mesh_to_spawn);
    ClassDB::bind_method(D_METHOD("get_mesh_to_spawn"), &UnderGenMeshSpawnerNode::get_mesh_to_spawn);

    ClassDB::bind_method(D_METHOD("set_material_override", "material"), &UnderGenMeshSpawnerNode::set_material_override);
    ClassDB::bind_method(D_METHOD("get_material_override"), &UnderGenMeshSpawnerNode::get_material_override);

    ClassDB::bind_method(D_METHOD("set_spawn_probability", "prob"), &UnderGenMeshSpawnerNode::set_spawn_probability);
    ClassDB::bind_method(D_METHOD("get_spawn_probability"), &UnderGenMeshSpawnerNode::get_spawn_probability);

    ClassDB::bind_method(D_METHOD("set_instances_per_point", "instances"), &UnderGenMeshSpawnerNode::set_instances_per_point);
    ClassDB::bind_method(D_METHOD("get_instances_per_point"), &UnderGenMeshSpawnerNode::get_instances_per_point);

    ClassDB::bind_method(D_METHOD("set_max_distance_from_point", "dist"), &UnderGenMeshSpawnerNode::set_max_distance_from_point);
    ClassDB::bind_method(D_METHOD("get_max_distance_from_point"), &UnderGenMeshSpawnerNode::get_max_distance_from_point);

    ClassDB::bind_method(D_METHOD("set_min_y_offset", "offset"), &UnderGenMeshSpawnerNode::set_min_y_offset);
    ClassDB::bind_method(D_METHOD("get_min_y_offset"), &UnderGenMeshSpawnerNode::get_min_y_offset);

    ClassDB::bind_method(D_METHOD("set_max_y_offset", "offset"), &UnderGenMeshSpawnerNode::set_max_y_offset);
    ClassDB::bind_method(D_METHOD("get_max_y_offset"), &UnderGenMeshSpawnerNode::get_max_y_offset);

    ClassDB::bind_method(D_METHOD("set_min_rotation", "rot"), &UnderGenMeshSpawnerNode::set_min_rotation);
    ClassDB::bind_method(D_METHOD("get_min_rotation"), &UnderGenMeshSpawnerNode::get_min_rotation);

    ClassDB::bind_method(D_METHOD("set_max_rotation", "rot"), &UnderGenMeshSpawnerNode::set_max_rotation);
    ClassDB::bind_method(D_METHOD("get_max_rotation"), &UnderGenMeshSpawnerNode::get_max_rotation);

    ClassDB::bind_method(D_METHOD("set_min_scale", "scale"), &UnderGenMeshSpawnerNode::set_min_scale);
    ClassDB::bind_method(D_METHOD("get_min_scale"), &UnderGenMeshSpawnerNode::get_min_scale);

    ClassDB::bind_method(D_METHOD("set_max_scale", "scale"), &UnderGenMeshSpawnerNode::set_max_scale);
    ClassDB::bind_method(D_METHOD("get_max_scale"), &UnderGenMeshSpawnerNode::get_max_scale);

    ClassDB::bind_method(D_METHOD("set_random_seed", "seed"), &UnderGenMeshSpawnerNode::set_random_seed);
    ClassDB::bind_method(D_METHOD("get_random_seed"), &UnderGenMeshSpawnerNode::get_random_seed);

    ClassDB::bind_method(D_METHOD("set_consume_points", "consume"), &UnderGenMeshSpawnerNode::set_consume_points);
    ClassDB::bind_method(D_METHOD("get_consume_points"), &UnderGenMeshSpawnerNode::get_consume_points);

    ClassDB::bind_method(D_METHOD("set_align_with_normal", "align"), &UnderGenMeshSpawnerNode::set_align_with_normal);
    ClassDB::bind_method(D_METHOD("get_align_with_normal"), &UnderGenMeshSpawnerNode::get_align_with_normal);

    ClassDB::bind_method(D_METHOD("set_cast_shadows", "shadows"), &UnderGenMeshSpawnerNode::set_cast_shadows);
    ClassDB::bind_method(D_METHOD("get_cast_shadows"), &UnderGenMeshSpawnerNode::get_cast_shadows);

    ClassDB::bind_method(D_METHOD("set_chunk_size", "size"), &UnderGenMeshSpawnerNode::set_chunk_size);
    ClassDB::bind_method(D_METHOD("get_chunk_size"), &UnderGenMeshSpawnerNode::get_chunk_size);

    ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "mesh_to_spawn", PROPERTY_HINT_RESOURCE_TYPE, "Mesh"), "set_mesh_to_spawn", "get_mesh_to_spawn");
    ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "material_override", PROPERTY_HINT_RESOURCE_TYPE, "Material"), "set_material_override", "get_material_override");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "spawn_probability", PROPERTY_HINT_RANGE, "0,1,0.01"), "set_spawn_probability", "get_spawn_probability");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "instances_per_point"), "set_instances_per_point", "get_instances_per_point");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "max_distance_from_point"), "set_max_distance_from_point", "get_max_distance_from_point");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "min_y_offset"), "set_min_y_offset", "get_min_y_offset");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "max_y_offset"), "set_max_y_offset", "get_max_y_offset");
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3, "min_rotation"), "set_min_rotation", "get_min_rotation");
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3, "max_rotation"), "set_max_rotation", "get_max_rotation");
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3, "min_scale"), "set_min_scale", "get_min_scale");
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3, "max_scale"), "set_max_scale", "get_max_scale");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "random_seed"), "set_random_seed", "get_random_seed");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "consume_points"), "set_consume_points", "get_consume_points");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "align_with_normal"), "set_align_with_normal", "get_align_with_normal");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "cast_shadows", PROPERTY_HINT_ENUM, "Off,On,Double-Sided,Shadows Only"), "set_cast_shadows", "get_cast_shadows");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "chunk_size"), "set_chunk_size", "get_chunk_size");
}

void UnderGenMeshSpawnerNode::set_mesh_to_spawn(const Ref<Mesh> &p_mesh) { mesh_to_spawn = p_mesh; }
Ref<Mesh> UnderGenMeshSpawnerNode::get_mesh_to_spawn() const { return mesh_to_spawn; }

void UnderGenMeshSpawnerNode::set_material_override(const Ref<Material> &p_material) { material_override = p_material; }
Ref<Material> UnderGenMeshSpawnerNode::get_material_override() const { return material_override; }

void UnderGenMeshSpawnerNode::set_spawn_probability(float p_prob) { spawn_probability = Math::clamp(p_prob, 0.0f, 1.0f); }
float UnderGenMeshSpawnerNode::get_spawn_probability() const { return spawn_probability; }

void UnderGenMeshSpawnerNode::set_instances_per_point(int p_instances) { instances_per_point = p_instances < 1 ? 1 : p_instances; }
int UnderGenMeshSpawnerNode::get_instances_per_point() const { return instances_per_point; }

void UnderGenMeshSpawnerNode::set_max_distance_from_point(float p_dist) { max_distance_from_point = p_dist < 0.0f ? 0.0f : p_dist; }
float UnderGenMeshSpawnerNode::get_max_distance_from_point() const { return max_distance_from_point; }

void UnderGenMeshSpawnerNode::set_min_y_offset(float p_offset) { min_y_offset = p_offset; }
float UnderGenMeshSpawnerNode::get_min_y_offset() const { return min_y_offset; }

void UnderGenMeshSpawnerNode::set_max_y_offset(float p_offset) { max_y_offset = p_offset; }
float UnderGenMeshSpawnerNode::get_max_y_offset() const { return max_y_offset; }

void UnderGenMeshSpawnerNode::set_min_rotation(const Vector3 &p_rot) { min_rotation = p_rot; }
Vector3 UnderGenMeshSpawnerNode::get_min_rotation() const { return min_rotation; }

void UnderGenMeshSpawnerNode::set_max_rotation(const Vector3 &p_rot) { max_rotation = p_rot; }
Vector3 UnderGenMeshSpawnerNode::get_max_rotation() const { return max_rotation; }

void UnderGenMeshSpawnerNode::set_min_scale(const Vector3 &p_scale) { min_scale = p_scale; }
Vector3 UnderGenMeshSpawnerNode::get_min_scale() const { return min_scale; }

void UnderGenMeshSpawnerNode::set_max_scale(const Vector3 &p_scale) { max_scale = p_scale; }
Vector3 UnderGenMeshSpawnerNode::get_max_scale() const { return max_scale; }

void UnderGenMeshSpawnerNode::set_random_seed(int64_t p_seed) { random_seed = p_seed; rng->set_seed(p_seed); }
int64_t UnderGenMeshSpawnerNode::get_random_seed() const { return random_seed; }

void UnderGenMeshSpawnerNode::set_consume_points(bool p_consume) { consume_points = p_consume; }
bool UnderGenMeshSpawnerNode::get_consume_points() const { return consume_points; }

void UnderGenMeshSpawnerNode::set_align_with_normal(bool p_align) { align_with_normal = p_align; }
bool UnderGenMeshSpawnerNode::get_align_with_normal() const { return align_with_normal; }

void UnderGenMeshSpawnerNode::set_cast_shadows(int p_shadows) { cast_shadows = p_shadows; }
int UnderGenMeshSpawnerNode::get_cast_shadows() const { return cast_shadows; }

void UnderGenMeshSpawnerNode::set_chunk_size(float p_size) { chunk_size = p_size; }
float UnderGenMeshSpawnerNode::get_chunk_size() const { return chunk_size; }

void UnderGenMeshSpawnerNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    outputs[0] = inputs.get(0, Ref<UnderGenPointSet>());
}

void UnderGenMeshSpawnerNode::execute_with_parent(const Dictionary &inputs, Dictionary &outputs, Node *target_parent) {
    if (!target_parent) {
        UtilityFunctions::printerr("UnderGenMeshSpawnerNode: No target parent node.");
        return;
    }

    Ref<UnderGenPointSet> point_set = inputs.get(0, Ref<UnderGenPointSet>());
    if (point_set.is_null()) {
        UtilityFunctions::printerr("UnderGenMeshSpawnerNode: No PointSet on port 0.");
        return;
    }

    if (mesh_to_spawn.is_null()) {
        UtilityFunctions::printerr("UnderGenMeshSpawnerNode: No mesh assigned.");
        return;
    }

    rng->set_seed(random_seed);

    const auto& raw = point_set->get_points_raw();
    if (raw.empty()) {
        outputs[0] = point_set;
        return;
    }

    std::map<std::tuple<int, int, int>, std::vector<Transform3D>> chunk_map;
    std::vector<size_t> spawned_indices;

    for (size_t i = 0; i < raw.size(); ++i) {
        const auto& p = raw[i];

        // Probabilistic rejection
        if (rng->randf() > spawn_probability * p.density) continue;

        bool spawned_any = false;

        for (int j = 0; j < instances_per_point; ++j) {
            Transform3D xform = p.transform;

            // 1. Align with normal if requested
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
            }

            // 2. Compute local offsets
            Vector3 offset;
            if (max_distance_from_point > 0.0f) {
                float theta = rng->randf() * Math_TAU;
                float r = Math::sqrt(rng->randf()) * max_distance_from_point;
                offset.x = Math::cos(theta) * r;
                offset.z = Math::sin(theta) * r;
            }
            if (min_y_offset != max_y_offset) {
                offset.y = rng->randf() * (max_y_offset - min_y_offset) + min_y_offset;
            } else {
                offset.y = min_y_offset;
            }

            // Translate transform by offset in local coordinates of the aligned basis
            xform.origin += xform.basis.xform(offset);

            // 3. Rotation variation
            Vector3 rot;
            rot.x = Math::deg_to_rad(rng->randf() * (max_rotation.x - min_rotation.x) + min_rotation.x);
            rot.y = Math::deg_to_rad(rng->randf() * (max_rotation.y - min_rotation.y) + min_rotation.y);
            rot.z = Math::deg_to_rad(rng->randf() * (max_rotation.z - min_rotation.z) + min_rotation.z);

            Basis rot_basis = Basis::from_euler(rot);
            xform.basis = xform.basis * rot_basis;

            // 4. Scale variation
            Vector3 scale;
            scale.x = rng->randf() * (max_scale.x - min_scale.x) + min_scale.x;
            scale.y = rng->randf() * (max_scale.y - min_scale.y) + min_scale.y;
            scale.z = rng->randf() * (max_scale.z - min_scale.z) + min_scale.z;

            xform.basis = xform.basis.scaled(scale);

            // Determine chunk coordinate
            int cx = 0, cy = 0, cz = 0;
            if (chunk_size > 0.0f) {
                cx = (int)Math::floor(xform.origin.x / chunk_size);
                cy = (int)Math::floor(xform.origin.y / chunk_size);
                cz = (int)Math::floor(xform.origin.z / chunk_size);
            }
            chunk_map[std::make_tuple(cx, cy, cz)].push_back(xform);
            spawned_any = true;
        }

        if (spawned_any) {
            spawned_indices.push_back(i);
        }
    }

    size_t total_spawned = 0;

    for (const auto& pair : chunk_map) {
        const auto& chunk_coords = pair.first;
        const auto& transforms = pair.second;
        if (transforms.empty()) continue;

        Ref<MultiMesh> mm;
        mm.instantiate();
        mm->set_transform_format(MultiMesh::TRANSFORM_3D);
        mm->set_mesh(mesh_to_spawn);
        mm->set_instance_count((int)transforms.size());

        for (int i = 0; i < (int)transforms.size(); ++i) {
            mm->set_instance_transform(i, transforms[i]);
        }

        MultiMeshInstance3D *mmi = memnew(MultiMeshInstance3D);
        mmi->set_multimesh(mm);
        
        int cx = std::get<0>(chunk_coords);
        int cy = std::get<1>(chunk_coords);
        int cz = std::get<2>(chunk_coords);

        if (chunk_size > 0.0f) {
            mmi->set_name("MeshDecoration_" + get_name() + "_" + String::num_int64(cx) + "_" + String::num_int64(cy) + "_" + String::num_int64(cz));
        } else {
            mmi->set_name("MeshDecoration_" + get_name());
        }

        if (material_override.is_valid()) {
            mmi->set_material_override(material_override);
        }
        mmi->set_cast_shadows_setting((GeometryInstance3D::ShadowCastingSetting)cast_shadows);

        target_parent->add_child(mmi);
        total_spawned += transforms.size();
    }

    if (consume_points && !spawned_indices.empty()) {
        std::sort(spawned_indices.begin(), spawned_indices.end(), std::greater<size_t>());
        for (size_t idx : spawned_indices) {
            point_set->remove_point((int)idx);
        }
    }

    UtilityFunctions::print("UnderGenMeshSpawnerNode: Spawned ", chunk_map.size(), " MultiMesh chunks, total ", total_spawned, " instances.");
    outputs[0] = point_set;
}

} // namespace godot
