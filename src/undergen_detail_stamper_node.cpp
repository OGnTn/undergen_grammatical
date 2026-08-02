#include "undergen_detail_stamper_node.h"
#include "density_grid.h"
#include <godot_cpp/classes/file_access.hpp>
#include <godot_cpp/variant/utility_functions.hpp>
#include <godot_cpp/core/math.hpp>
#include <algorithm>
#include <set>

namespace godot {

UnderGenDetailStamperNode::UnderGenDetailStamperNode() {
    rng.instantiate();
    align_yaw_only = false;
}

UnderGenDetailStamperNode::~UnderGenDetailStamperNode() {
    _clear_vox_cache();
}

void UnderGenDetailStamperNode::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_vox_path", "path"), &UnderGenDetailStamperNode::set_vox_path);
    ClassDB::bind_method(D_METHOD("get_vox_path"), &UnderGenDetailStamperNode::get_vox_path);
    ClassDB::bind_method(D_METHOD("set_probability", "prob"), &UnderGenDetailStamperNode::set_probability);
    ClassDB::bind_method(D_METHOD("get_probability"), &UnderGenDetailStamperNode::get_probability);
    ClassDB::bind_method(D_METHOD("set_align_with_normal", "align"), &UnderGenDetailStamperNode::set_align_with_normal);
    ClassDB::bind_method(D_METHOD("get_align_with_normal"), &UnderGenDetailStamperNode::get_align_with_normal);
    ClassDB::bind_method(D_METHOD("set_align_yaw_only", "enabled"), &UnderGenDetailStamperNode::set_align_yaw_only);
    ClassDB::bind_method(D_METHOD("get_align_yaw_only"), &UnderGenDetailStamperNode::get_align_yaw_only);
    ClassDB::bind_method(D_METHOD("set_random_y_rotation", "rand"), &UnderGenDetailStamperNode::set_random_y_rotation);
    ClassDB::bind_method(D_METHOD("get_random_y_rotation"), &UnderGenDetailStamperNode::get_random_y_rotation);
    ClassDB::bind_method(D_METHOD("set_vox_inverse_density", "enabled"), &UnderGenDetailStamperNode::set_vox_inverse_density);
    ClassDB::bind_method(D_METHOD("get_vox_inverse_density"), &UnderGenDetailStamperNode::get_vox_inverse_density);
    ClassDB::bind_method(D_METHOD("set_pivot_preset", "preset"), &UnderGenDetailStamperNode::set_pivot_preset);
    ClassDB::bind_method(D_METHOD("get_pivot_preset"), &UnderGenDetailStamperNode::get_pivot_preset);
    ClassDB::bind_method(D_METHOD("set_custom_pivot", "pivot"), &UnderGenDetailStamperNode::set_custom_pivot);
    ClassDB::bind_method(D_METHOD("get_custom_pivot"), &UnderGenDetailStamperNode::get_custom_pivot);
    ClassDB::bind_method(D_METHOD("set_random_seed", "seed"), &UnderGenDetailStamperNode::set_random_seed);
    ClassDB::bind_method(D_METHOD("get_random_seed"), &UnderGenDetailStamperNode::get_random_seed);

    ClassDB::bind_method(D_METHOD("set_vox_spawn_entries", "entries"), &UnderGenDetailStamperNode::set_vox_spawn_entries);
    ClassDB::bind_method(D_METHOD("get_vox_spawn_entries"), &UnderGenDetailStamperNode::get_vox_spawn_entries);
    ClassDB::bind_method(D_METHOD("set_vox_material_entries", "entries"), &UnderGenDetailStamperNode::set_vox_material_entries);
    ClassDB::bind_method(D_METHOD("get_vox_material_entries"), &UnderGenDetailStamperNode::get_vox_material_entries);

    ClassDB::bind_method(D_METHOD("_on_entries_changed"), &UnderGenDetailStamperNode::_on_entries_changed);

    BIND_ENUM_CONSTANT(PIVOT_BOTTOM_CENTER);
    BIND_ENUM_CONSTANT(PIVOT_TOP_CENTER);
    BIND_ENUM_CONSTANT(PIVOT_CENTER);
    BIND_ENUM_CONSTANT(PIVOT_CUSTOM);

    ADD_PROPERTY(PropertyInfo(Variant::STRING, "vox_path", PROPERTY_HINT_FILE, "*.vox"), "set_vox_path", "get_vox_path");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "probability", PROPERTY_HINT_RANGE, "0.0,1.0,0.05"), "set_probability", "get_probability");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "align_with_normal"), "set_align_with_normal", "get_align_with_normal");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "align_yaw_only"), "set_align_yaw_only", "get_align_yaw_only");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "random_y_rotation"), "set_random_y_rotation", "get_random_y_rotation");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "vox_inverse_density"), "set_vox_inverse_density", "get_vox_inverse_density");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "pivot_preset", PROPERTY_HINT_ENUM, "Bottom Center,Top Center,Center,Custom"), "set_pivot_preset", "get_pivot_preset");
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3, "custom_pivot"), "set_custom_pivot", "get_custom_pivot");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "random_seed"), "set_random_seed", "get_random_seed");

    ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "vox_spawn_entries", PROPERTY_HINT_TYPE_STRING,
                              vformat("%d/%d:%s", Variant::OBJECT, Variant::NIL, "VoxSpawnEntry")),
                 "set_vox_spawn_entries", "get_vox_spawn_entries");
    ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "vox_material_entries", PROPERTY_HINT_TYPE_STRING,
                              vformat("%d/%d:%s", Variant::OBJECT, Variant::NIL, "VoxMaterialEntry")),
                 "set_vox_material_entries", "get_vox_material_entries");
}

void UnderGenDetailStamperNode::set_vox_path(const String &p_path) { vox_path = p_path; emit_changed(); }
String UnderGenDetailStamperNode::get_vox_path() const { return vox_path; }

void UnderGenDetailStamperNode::set_probability(float p_prob) { probability = p_prob; emit_changed(); }
float UnderGenDetailStamperNode::get_probability() const { return probability; }

void UnderGenDetailStamperNode::set_align_with_normal(bool p_align) { align_with_normal = p_align; emit_changed(); }
bool UnderGenDetailStamperNode::get_align_with_normal() const { return align_with_normal; }

void UnderGenDetailStamperNode::set_align_yaw_only(bool p_enabled) { align_yaw_only = p_enabled; emit_changed(); }
bool UnderGenDetailStamperNode::get_align_yaw_only() const { return align_yaw_only; }

void UnderGenDetailStamperNode::set_random_y_rotation(bool p_rand) { random_y_rotation = p_rand; emit_changed(); }
bool UnderGenDetailStamperNode::get_random_y_rotation() const { return random_y_rotation; }

void UnderGenDetailStamperNode::set_vox_inverse_density(bool p_enabled) { vox_inverse_density = p_enabled; emit_changed(); }
bool UnderGenDetailStamperNode::get_vox_inverse_density() const { return vox_inverse_density; }

void UnderGenDetailStamperNode::set_pivot_preset(PivotPreset p_preset) { pivot_preset = p_preset; emit_changed(); }
UnderGenDetailStamperNode::PivotPreset UnderGenDetailStamperNode::get_pivot_preset() const { return pivot_preset; }

void UnderGenDetailStamperNode::set_custom_pivot(const Vector3 &p_pivot) { custom_pivot = p_pivot; emit_changed(); }
Vector3 UnderGenDetailStamperNode::get_custom_pivot() const { return custom_pivot; }

void UnderGenDetailStamperNode::set_random_seed(int p_seed) { random_seed = p_seed; emit_changed(); }
int UnderGenDetailStamperNode::get_random_seed() const { return random_seed; }

void UnderGenDetailStamperNode::set_vox_spawn_entries(const TypedArray<VoxSpawnEntry> &p_entries) {
    for (int i = 0; i < vox_spawn_entries.size(); ++i) {
        Ref<VoxSpawnEntry> entry = vox_spawn_entries[i];
        if (entry.is_valid() && entry->is_connected("changed", Callable(this, "_on_entries_changed"))) {
            entry->disconnect("changed", Callable(this, "_on_entries_changed"));
        }
    }
    vox_spawn_entries = p_entries;
    for (int i = 0; i < vox_spawn_entries.size(); ++i) {
        Ref<VoxSpawnEntry> entry = vox_spawn_entries[i];
        if (entry.is_valid() && !entry->is_connected("changed", Callable(this, "_on_entries_changed"))) {
            entry->connect("changed", Callable(this, "_on_entries_changed"));
        }
    }
    _rebuild_maps_from_entries();
    emit_changed();
}

TypedArray<VoxSpawnEntry> UnderGenDetailStamperNode::get_vox_spawn_entries() const { return vox_spawn_entries; }

void UnderGenDetailStamperNode::set_vox_material_entries(const TypedArray<VoxMaterialEntry> &p_entries) {
    for (int i = 0; i < vox_material_entries.size(); ++i) {
        Ref<VoxMaterialEntry> entry = vox_material_entries[i];
        if (entry.is_valid() && entry->is_connected("changed", Callable(this, "_on_entries_changed"))) {
            entry->disconnect("changed", Callable(this, "_on_entries_changed"));
        }
    }
    vox_material_entries = p_entries;
    for (int i = 0; i < vox_material_entries.size(); ++i) {
        Ref<VoxMaterialEntry> entry = vox_material_entries[i];
        if (entry.is_valid() && !entry->is_connected("changed", Callable(this, "_on_entries_changed"))) {
            entry->connect("changed", Callable(this, "_on_entries_changed"));
        }
    }
    _rebuild_maps_from_entries();
    emit_changed();
}

TypedArray<VoxMaterialEntry> UnderGenDetailStamperNode::get_vox_material_entries() const { return vox_material_entries; }

void UnderGenDetailStamperNode::_rebuild_maps_from_entries() {
    vox_spawn_map.clear();
    for (int i = 0; i < vox_spawn_entries.size(); ++i) {
        Ref<VoxSpawnEntry> entry = vox_spawn_entries[i];
        if (entry.is_valid()) {
            vox_spawn_map[entry->get_palette_index()] = entry->get_spawn_type();
        }
    }
    vox_material_map.clear();
    for (int i = 0; i < vox_material_entries.size(); ++i) {
        Ref<VoxMaterialEntry> entry = vox_material_entries[i];
        if (entry.is_valid()) {
            vox_material_map[entry->get_palette_index()] = entry->get_material_id();
        }
    }
}

void UnderGenDetailStamperNode::_on_entries_changed() {
    _rebuild_maps_from_entries();
    emit_changed();
}

const ogt_vox_scene* UnderGenDetailStamperNode::_load_vox(const String &path) {
    if (path.is_empty()) return nullptr;
    auto it = vox_cache.find(path);
    if (it != vox_cache.end()) return it->second;

    Ref<FileAccess> file = FileAccess::open(path, FileAccess::READ);
    if (file.is_null()) {
        UtilityFunctions::printerr("UnderGenDetailStamperNode: Failed to open vox: ", path);
        return nullptr;
    }
    uint64_t len = file->get_length();
    PackedByteArray buf = file->get_buffer(len);
    const ogt_vox_scene* scene = ogt_vox_read_scene(buf.ptr(), (uint32_t)len);
    if (!scene) {
        UtilityFunctions::printerr("UnderGenDetailStamperNode: Failed to parse vox: ", path);
        return nullptr;
    }
    vox_cache[path] = scene;
    return scene;
}

void UnderGenDetailStamperNode::_clear_vox_cache() {
    for (auto const& [path, scene] : vox_cache) {
        ogt_vox_destroy_scene(scene);
    }
    vox_cache.clear();
}

void UnderGenDetailStamperNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    // Port 0: context
    // Port 1: PointSet
    Dictionary context = inputs.get(0, Dictionary());
    Ref<UnderGenPointSet> point_set = inputs.get(1, Ref<UnderGenPointSet>());

    if (context.is_empty() || point_set.is_null()) {
        outputs[0] = context;
        return;
    }

    Ref<DensityGrid> grid = context.get("grid", Ref<DensityGrid>());
    if (grid.is_null() || vox_path.is_empty()) {
        outputs[0] = context;
        return;
    }

    const ogt_vox_scene* scene = _load_vox(vox_path);
    if (!scene) {
        outputs[0] = context;
        return;
    }

    // 1. Gather visible model dimensions & bbox to calculate pivot
    int scene_min_x = INT_MAX, scene_max_x = INT_MIN;
    int scene_min_y = INT_MAX, scene_max_y = INT_MIN;
    int scene_min_z = INT_MAX, scene_max_z = INT_MIN;
    bool has_visible = false;

    for (uint32_t i = 0; i < scene->num_instances; ++i) {
        const ogt_vox_instance& inst = scene->instances[i];
        if (inst.hidden) continue;
        const ogt_vox_model* model = scene->models[inst.model_index];
        if (!model) continue;

        has_visible = true;

        int offset_x = (int)inst.transform.m30;
        int offset_y = (int)inst.transform.m31;
        int offset_z = (int)inst.transform.m32;

        scene_min_x = std::min(scene_min_x, offset_x);
        scene_max_x = std::max(scene_max_x, offset_x + (int)model->size_x);

        scene_min_y = std::min(scene_min_y, offset_z);
        scene_max_y = std::max(scene_max_y, offset_z + (int)model->size_z);

        scene_min_z = std::min(scene_min_z, offset_y);
        scene_max_z = std::max(scene_max_z, offset_y + (int)model->size_y);
    }

    if (!has_visible) {
        outputs[0] = context;
        return;
    }

    int scene_size_x = scene_max_x - scene_min_x;
    int scene_size_y = scene_max_y - scene_min_y;
    int scene_size_z = scene_max_z - scene_min_z;

    Vector3 pivot;
    switch (pivot_preset) {
        case PIVOT_BOTTOM_CENTER:
            pivot = Vector3(scene_min_x + scene_size_x / 2.0f, scene_min_y, scene_min_z + scene_size_z / 2.0f);
            break;
        case PIVOT_TOP_CENTER:
            pivot = Vector3(scene_min_x + scene_size_x / 2.0f, scene_min_y + scene_size_y, scene_min_z + scene_size_z / 2.0f);
            break;
        case PIVOT_CENTER:
            pivot = Vector3(scene_min_x + scene_size_x / 2.0f, scene_min_y + scene_size_y / 2.0f, scene_min_z + scene_size_z / 2.0f);
            break;
        case PIVOT_CUSTOM:
            pivot = Vector3(scene_min_x + custom_pivot.x, scene_min_y + custom_pivot.y, scene_min_z + custom_pivot.z);
            break;
    }

    int gsx = grid->get_grid_size_x();
    int gsy = grid->get_grid_size_y();
    int gsz = grid->get_grid_size_z();

    constexpr float OPEN = 0.0f;
    constexpr float SOLID = 1.0f;

    // Load active voxel size
    float v_size = context.get("voxel_size", 1.0f);

    rng->set_seed(random_seed);

    // Rebuild entries list to ensure connection/palette values are correct
    for (int i = 0; i < vox_spawn_entries.size(); ++i) {
        Ref<VoxSpawnEntry> entry = vox_spawn_entries[i];
        if (entry.is_valid() && !entry->is_connected("changed", Callable(this, "_on_entries_changed"))) {
            entry->connect("changed", Callable(this, "_on_entries_changed"));
        }
    }
    for (int i = 0; i < vox_material_entries.size(); ++i) {
        Ref<VoxMaterialEntry> entry = vox_material_entries[i];
        if (entry.is_valid() && !entry->is_connected("changed", Callable(this, "_on_entries_changed"))) {
            entry->connect("changed", Callable(this, "_on_entries_changed"));
        }
    }
    _rebuild_maps_from_entries();

    // Get current spawn array to append to
    Array spawns_array = context.get("vox_spawns", Array());

    const auto& raw_points = point_set->get_points_raw();
    int stamp_count = 0;

    for (const auto& p : raw_points) {
        if (rng->randf() > probability * p.density) continue;

        // Voxel coordinate from point's world position
        Vector3i grid_origin(
            (int)Math::round(p.transform.origin.x / v_size),
            (int)Math::round(p.transform.origin.y / v_size),
            (int)Math::round(p.transform.origin.z / v_size)
        );

        // Rotation basis
        Basis rotation_basis;
        if (align_with_normal) {
            Vector3 normal = p.attributes.get("normal", Vector3(0, 1, 0));
            if (align_yaw_only) {
                Vector3 flat_normal(normal.x, 0.0f, normal.z);
                if (flat_normal.length_squared() > 0.001f) {
                    flat_normal.normalize();
                    float yaw = Math::atan2(flat_normal.x, flat_normal.z);
                    rotation_basis = Basis(Vector3(0, 1, 0), yaw);
                }
            } else {
                Vector3 up = Vector3(0, 1, 0);
                if (normal.is_equal_approx(-up)) {
                    rotation_basis = Basis(Vector3(1, 0, 0), Math_PI);
                } else if (!normal.is_equal_approx(up)) {
                    Vector3 axis = up.cross(normal).normalized();
                    float angle = Math::acos(up.dot(normal));
                    rotation_basis = Basis(axis, angle);
                }
            }

            if (random_y_rotation) {
                float angle = rng->randf() * Math_TAU;
                Vector3 rot_axis = align_yaw_only ? Vector3(0, 1, 0) : normal;
                rotation_basis = rotation_basis.rotated(rot_axis, angle);
            }
        } else {
            if (random_y_rotation) {
                float angle = rng->randf() * Math_TAU;
                rotation_basis = rotation_basis.rotated(Vector3(0, 1, 0), angle);
            }
        }

        // Loop over instances & stamp voxels
        for (uint32_t idx = 0; idx < scene->num_instances; ++idx) {
            const ogt_vox_instance& inst = scene->instances[idx];
            if (inst.hidden) continue;
            const ogt_vox_model* model = scene->models[inst.model_index];
            if (!model) continue;

            int offset_x = (int)inst.transform.m30;
            int offset_y = (int)inst.transform.m31;
            int offset_z = (int)inst.transform.m32;

            for (int vx = 0; vx < (int)model->size_x; ++vx) {
                for (int vy = 0; vy < (int)model->size_y; ++vy) {
                    for (int vz = 0; vz < (int)model->size_z; ++vz) {
                        uint32_t vi = vx + (vy * model->size_x) + (vz * model->size_x * model->size_y);
                        uint8_t ci = model->voxel_data[vi];
                        if (ci == 0) continue; // Skip air

                        // Local coordinate in MV -> Godot space
                        int rx = offset_x + vx;
                        int ry = offset_z + vz;
                        int rz = offset_y + vy;

                        Vector3 local_pos = Vector3(rx, ry, rz) - pivot;
                        Vector3 rotated_pos = rotation_basis.xform(local_pos);

                        Vector3i gp = grid_origin + Vector3i(
                            (int)Math::round(rotated_pos.x),
                            (int)Math::round(rotated_pos.y),
                            (int)Math::round(rotated_pos.z)
                        );

                        if (gp.x < 0 || gp.x >= gsx || gp.y < 0 || gp.y >= gsy || gp.z < 0 || gp.z >= gsz) {
                            continue;
                        }

                        if (vox_inverse_density) {
                            grid->set_cell(gp, OPEN);
                        } else {
                            Variant key_int = (int)ci;
                            Variant key_str = String::num_int64(ci);

                            if (vox_spawn_map.has(key_int)) {
                                grid->set_cell(gp, OPEN);
                                Dictionary spawn_d;
                                spawn_d["position"] = Vector3(gp.x, gp.y, gp.z);
                                spawn_d["palette_index"] = (int)ci;
                                spawn_d["spawn_type"] = vox_spawn_map[key_int];
                                spawns_array.append(spawn_d);
                            } else if (vox_spawn_map.has(key_str)) {
                                grid->set_cell(gp, OPEN);
                                Dictionary spawn_d;
                                spawn_d["position"] = Vector3(gp.x, gp.y, gp.z);
                                spawn_d["palette_index"] = (int)ci;
                                spawn_d["spawn_type"] = vox_spawn_map[key_str];
                                spawns_array.append(spawn_d);
                            } else {
                                grid->set_cell(gp, SOLID);
                                if (vox_material_map.has(key_int)) {
                                    grid->set_material_id(gp, (int)vox_material_map[key_int]);
                                } else if (vox_material_map.has(key_str)) {
                                    grid->set_material_id(gp, (int)vox_material_map[key_str]);
                                }
                            }
                        }
                    }
                }
            }
        }
        stamp_count++;
    }

    context["vox_spawns"] = spawns_array;

    // Output material thicknesses and stepped properties to the generation context
    Dictionary thicknesses = context.get("material_thicknesses", Dictionary());
    Dictionary stepped_dict = context.get("material_stepped", Dictionary());
    for (int i = 0; i < vox_material_entries.size(); ++i) {
        Ref<VoxMaterialEntry> entry = vox_material_entries[i];
        if (entry.is_valid()) {
            int mat_id = entry->get_material_id();
            if (mat_id > 0) {
                thicknesses[mat_id] = entry->get_thickness();
                stepped_dict[mat_id] = entry->is_stepped();
            }
        }
    }
    context["material_thicknesses"] = thicknesses;
    context["material_stepped"] = stepped_dict;

    UtilityFunctions::print("UnderGenDetailStamperNode: Stamped details at ", stamp_count, " / ", raw_points.size(), " points using '", vox_path, "'.");

    outputs[0] = context;
}

} // namespace godot
