#include "undergen_vox_stamp_node.h"
#include <godot_cpp/classes/file_access.hpp>
#include <godot_cpp/variant/utility_functions.hpp>
#include <godot_cpp/core/math.hpp>
#include <set>

#define OGT_VOX_IMPLEMENTATION
#include "ogt_vox.h"

namespace godot {

UnderGenVoxStampNode::UnderGenVoxStampNode() {}
UnderGenVoxStampNode::~UnderGenVoxStampNode() {
    _clear_vox_cache();
}

void UnderGenVoxStampNode::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_vox_spawn_map", "map"), &UnderGenVoxStampNode::set_vox_spawn_map);
    ClassDB::bind_method(D_METHOD("get_vox_spawn_map"), &UnderGenVoxStampNode::get_vox_spawn_map);
    ClassDB::bind_method(D_METHOD("set_vox_material_map", "map"), &UnderGenVoxStampNode::set_vox_material_map);
    ClassDB::bind_method(D_METHOD("get_vox_material_map"), &UnderGenVoxStampNode::get_vox_material_map);
    ClassDB::bind_method(D_METHOD("set_vox_inverse_density", "enabled"), &UnderGenVoxStampNode::set_vox_inverse_density);
    ClassDB::bind_method(D_METHOD("get_vox_inverse_density"), &UnderGenVoxStampNode::get_vox_inverse_density);
    ClassDB::bind_method(D_METHOD("set_connection_palette_index", "palette_index"), &UnderGenVoxStampNode::set_connection_palette_index);
    ClassDB::bind_method(D_METHOD("get_connection_palette_index"), &UnderGenVoxStampNode::get_connection_palette_index);
    ClassDB::bind_method(D_METHOD("set_exclude_from_smoothing", "exclude"), &UnderGenVoxStampNode::set_exclude_from_smoothing);
    ClassDB::bind_method(D_METHOD("get_exclude_from_smoothing"), &UnderGenVoxStampNode::get_exclude_from_smoothing);

    ClassDB::bind_method(D_METHOD("set_vox_spawn_entries", "entries"), &UnderGenVoxStampNode::set_vox_spawn_entries);
    ClassDB::bind_method(D_METHOD("get_vox_spawn_entries"), &UnderGenVoxStampNode::get_vox_spawn_entries);
    ClassDB::bind_method(D_METHOD("set_vox_material_entries", "entries"), &UnderGenVoxStampNode::set_vox_material_entries);
    ClassDB::bind_method(D_METHOD("get_vox_material_entries"), &UnderGenVoxStampNode::get_vox_material_entries);

    ClassDB::bind_method(D_METHOD("_on_entries_changed"), &UnderGenVoxStampNode::_on_entries_changed);

    // Keep hidden from editor but serialized
    ADD_PROPERTY(PropertyInfo(Variant::DICTIONARY, "vox_spawn_map", PROPERTY_HINT_NONE, "", PROPERTY_USAGE_NO_EDITOR), "set_vox_spawn_map", "get_vox_spawn_map");
    ADD_PROPERTY(PropertyInfo(Variant::DICTIONARY, "vox_material_map", PROPERTY_HINT_NONE, "", PROPERTY_USAGE_NO_EDITOR), "set_vox_material_map", "get_vox_material_map");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "vox_inverse_density"), "set_vox_inverse_density", "get_vox_inverse_density");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "connection_palette_index"), "set_connection_palette_index", "get_connection_palette_index");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "exclude_from_smoothing"), "set_exclude_from_smoothing", "get_exclude_from_smoothing");

    // Exposed to editor
    ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "vox_spawn_entries", PROPERTY_HINT_TYPE_STRING,
                              vformat("%d/%d:%s", Variant::OBJECT, Variant::NIL, "VoxSpawnEntry")),
                 "set_vox_spawn_entries", "get_vox_spawn_entries");
    ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "vox_material_entries", PROPERTY_HINT_TYPE_STRING,
                              vformat("%d/%d:%s", Variant::OBJECT, Variant::NIL, "VoxMaterialEntry")),
                 "set_vox_material_entries", "get_vox_material_entries");
}

void UnderGenVoxStampNode::set_vox_spawn_map(const Dictionary &p_map) {
    vox_spawn_map = p_map;
    
    // Disconnect old entries
    for (int i = 0; i < vox_spawn_entries.size(); ++i) {
        Ref<VoxSpawnEntry> entry = vox_spawn_entries[i];
        if (entry.is_valid() && entry->is_connected("changed", Callable(this, "_on_entries_changed"))) {
            entry->disconnect("changed", Callable(this, "_on_entries_changed"));
        }
    }
    
    vox_spawn_entries.clear();
    Array keys = p_map.keys();
    for (int i = 0; i < keys.size(); ++i) {
        Variant key = keys[i];
        Ref<VoxSpawnEntry> entry;
        entry.instantiate();
        
        int idx = 0;
        if (key.get_type() == Variant::INT) {
            idx = (int)key;
        } else if (key.get_type() == Variant::STRING) {
            idx = ((String)key).to_int();
        }
        
        entry->set_palette_index(idx);
        entry->set_spawn_type(p_map[key]);
        
        entry->connect("changed", Callable(this, "_on_entries_changed"));
        vox_spawn_entries.append(entry);
    }
    
    emit_changed();
}

Dictionary UnderGenVoxStampNode::get_vox_spawn_map() const {
    return vox_spawn_map;
}

void UnderGenVoxStampNode::set_vox_material_map(const Dictionary &p_map) {
    vox_material_map = p_map;
    
    // Disconnect old entries
    for (int i = 0; i < vox_material_entries.size(); ++i) {
        Ref<VoxMaterialEntry> entry = vox_material_entries[i];
        if (entry.is_valid() && entry->is_connected("changed", Callable(this, "_on_entries_changed"))) {
            entry->disconnect("changed", Callable(this, "_on_entries_changed"));
        }
    }
    
    vox_material_entries.clear();
    Array keys = p_map.keys();
    for (int i = 0; i < keys.size(); ++i) {
        Variant key = keys[i];
        Ref<VoxMaterialEntry> entry;
        entry.instantiate();
        
        int idx = 0;
        if (key.get_type() == Variant::INT) {
            idx = (int)key;
        } else if (key.get_type() == Variant::STRING) {
            idx = ((String)key).to_int();
        }
        
        entry->set_palette_index(idx);
        entry->set_material_id(p_map[key]);
        
        entry->connect("changed", Callable(this, "_on_entries_changed"));
        vox_material_entries.append(entry);
    }
    
    emit_changed();
}

Dictionary UnderGenVoxStampNode::get_vox_material_map() const {
    return vox_material_map;
}

void UnderGenVoxStampNode::set_vox_spawn_entries(const TypedArray<VoxSpawnEntry> &p_entries) {
    // Disconnect old entries
    for (int i = 0; i < vox_spawn_entries.size(); ++i) {
        Ref<VoxSpawnEntry> entry = vox_spawn_entries[i];
        if (entry.is_valid() && entry->is_connected("changed", Callable(this, "_on_entries_changed"))) {
            entry->disconnect("changed", Callable(this, "_on_entries_changed"));
        }
    }
    
    vox_spawn_entries = p_entries;
    
    // Connect new entries
    for (int i = 0; i < vox_spawn_entries.size(); ++i) {
        Ref<VoxSpawnEntry> entry = vox_spawn_entries[i];
        if (entry.is_valid()) {
            if (!entry->is_connected("changed", Callable(this, "_on_entries_changed"))) {
                entry->connect("changed", Callable(this, "_on_entries_changed"));
            }
        }
    }
    
    _rebuild_maps_from_entries();
    emit_changed();
}

TypedArray<VoxSpawnEntry> UnderGenVoxStampNode::get_vox_spawn_entries() const {
    return vox_spawn_entries;
}

void UnderGenVoxStampNode::set_vox_material_entries(const TypedArray<VoxMaterialEntry> &p_entries) {
    // Disconnect old entries
    for (int i = 0; i < vox_material_entries.size(); ++i) {
        Ref<VoxMaterialEntry> entry = vox_material_entries[i];
        if (entry.is_valid() && entry->is_connected("changed", Callable(this, "_on_entries_changed"))) {
            entry->disconnect("changed", Callable(this, "_on_entries_changed"));
        }
    }
    
    vox_material_entries = p_entries;
    
    // Connect new entries
    for (int i = 0; i < vox_material_entries.size(); ++i) {
        Ref<VoxMaterialEntry> entry = vox_material_entries[i];
        if (entry.is_valid()) {
            if (!entry->is_connected("changed", Callable(this, "_on_entries_changed"))) {
                entry->connect("changed", Callable(this, "_on_entries_changed"));
            }
        }
    }
    
    _rebuild_maps_from_entries();
    emit_changed();
}

TypedArray<VoxMaterialEntry> UnderGenVoxStampNode::get_vox_material_entries() const {
    return vox_material_entries;
}

void UnderGenVoxStampNode::_rebuild_maps_from_entries() {
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

void UnderGenVoxStampNode::_on_entries_changed() {
    _rebuild_maps_from_entries();
    emit_changed();
}

void UnderGenVoxStampNode::set_vox_inverse_density(bool p_enabled) { vox_inverse_density = p_enabled; emit_changed(); }
bool UnderGenVoxStampNode::get_vox_inverse_density() const { return vox_inverse_density; }
void UnderGenVoxStampNode::set_connection_palette_index(int p_index) { connection_palette_index = p_index; emit_changed(); }
int UnderGenVoxStampNode::get_connection_palette_index() const { return connection_palette_index; }
void UnderGenVoxStampNode::set_exclude_from_smoothing(bool p_exclude) { exclude_from_smoothing = p_exclude; emit_changed(); }
bool UnderGenVoxStampNode::get_exclude_from_smoothing() const { return exclude_from_smoothing; }

const ogt_vox_scene* UnderGenVoxStampNode::_load_vox(const String &path) {
    if (path.is_empty()) return nullptr;
    auto it = vox_cache.find(path);
    if (it != vox_cache.end()) return it->second;

    Ref<FileAccess> file = FileAccess::open(path, FileAccess::READ);
    if (file.is_null()) {
        UtilityFunctions::printerr("UnderGenVoxStampNode: Failed to open vox: ", path);
        return nullptr;
    }
    uint64_t len = file->get_length();
    PackedByteArray buf = file->get_buffer(len);
    const ogt_vox_scene* scene = ogt_vox_read_scene(buf.ptr(), (uint32_t)len);
    if (!scene) {
        UtilityFunctions::printerr("UnderGenVoxStampNode: Failed to parse vox: ", path);
        return nullptr;
    }
    vox_cache[path] = scene;
    return scene;
}

void UnderGenVoxStampNode::_clear_vox_cache() {
    for (auto const& [path, scene] : vox_cache) {
        ogt_vox_destroy_scene(scene);
    }
    vox_cache.clear();
}

void UnderGenVoxStampNode::_stamp_vox(DensityGrid* grid, const ResolvedRoom &room,
                                       const ogt_vox_scene* scene,
                                       std::vector<Dictionary> &out_spawns,
                                       Array &out_connections) {
    if (!scene || !grid) return;

    UtilityFunctions::print("UnderGenVoxStampNode: Stamping '", room.vox_path, "'. Active vox_spawn_map: ", vox_spawn_map, " Active vox_material_map: ", vox_material_map);

    int gsx = grid->get_grid_size_x();
    int gsy = grid->get_grid_size_y();
    int gsz = grid->get_grid_size_z();
    Vector3i room_origin = room.position;

    constexpr float OPEN = 0.0f;
    constexpr float SOLID = 1.0f;

    std::set<uint8_t> unique_indices;

    // 1. Scan all instances to find the bounding box of the scene in Godot space
    int scene_min_x = INT_MAX, scene_max_x = INT_MIN;
    int scene_min_y = INT_MAX, scene_max_y = INT_MIN;
    int scene_min_z = INT_MAX, scene_max_z = INT_MIN;
    bool has_visible_instances = false;

    for (uint32_t i = 0; i < scene->num_instances; ++i) {
        const ogt_vox_instance& inst = scene->instances[i];
        if (inst.hidden) continue;
        const ogt_vox_model* model = scene->models[inst.model_index];
        if (!model) continue;

        has_visible_instances = true;

        int offset_x = (int)inst.transform.m30;
        int offset_y = (int)inst.transform.m31;
        int offset_z = (int)inst.transform.m32;

        int inst_min_x = offset_x;
        int inst_max_x = offset_x + model->size_x;
        int inst_min_y = offset_z; // MV Z -> Godot Y (vertical)
        int inst_max_y = offset_z + model->size_z;
        int inst_min_z = offset_y; // MV Y -> Godot Z (depth)
        int inst_max_z = offset_y + model->size_y;

        if (inst_min_x < scene_min_x) scene_min_x = inst_min_x;
        if (inst_max_x > scene_max_x) scene_max_x = inst_max_x;
        if (inst_min_y < scene_min_y) scene_min_y = inst_min_y;
        if (inst_max_y > scene_max_y) scene_max_y = inst_max_y;
        if (inst_min_z < scene_min_z) scene_min_z = inst_min_z;
        if (inst_max_z > scene_max_z) scene_max_z = inst_max_z;
    }

    if (!has_visible_instances) {
        UtilityFunctions::print("UnderGenVoxStampNode: No visible instances in scene, skipping.");
        return;
    }

    int scene_size_x = scene_max_x - scene_min_x;
    int scene_size_y = scene_max_y - scene_min_y; // Godot Y size (height)
    int scene_size_z = scene_max_z - scene_min_z; // Godot Z size (depth)

    // 2. Stamp each instance relative to the calculated bounding box
    for (uint32_t i = 0; i < scene->num_instances; ++i) {
        const ogt_vox_instance& inst = scene->instances[i];
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

                    if (ci != 0) {
                        unique_indices.insert(ci);
                    }

                    int rx = offset_x + vx;
                    int ry = offset_z + vz;
                    int rz = offset_y + vy;

                    // fx, fz: Centered horizontally in the room leaf
                    // fy: Aligned to room bottom (floor)
                    int fx = room_origin.x + (room.size.x - scene_size_x) / 2 + (rx - scene_min_x);
                    int fy = room_origin.y + (ry - scene_min_y);
                    int fz = room_origin.z + (room.size.z - scene_size_z) / 2 + (rz - scene_min_z);

                    if (fx < 0 || fx >= gsx || fy < 0 || fy >= gsy || fz < 0 || fz >= gsz) {
                        continue;
                    }
                    Vector3i gp(fx, fy, fz);

                    if (ci == 0) { // Air
                        if (!vox_inverse_density) grid->set_cell(gp, OPEN);
                    } else { // Solid
                        Variant key_int = (int)ci;
                        Variant key_str = String::num_int64(ci);
                        if (ci == connection_palette_index) {
                            grid->set_cell(gp, OPEN);
                            out_connections.append(Vector3(fx, fy, fz));
                        } else if (vox_spawn_map.has(key_int)) {
                            grid->set_cell(gp, OPEN);
                            Dictionary spawn_d;
                            spawn_d["position"] = Vector3(fx, fy, fz);
                            spawn_d["palette_index"] = (int)ci;
                            spawn_d["spawn_type"] = vox_spawn_map[key_int];
                            out_spawns.push_back(spawn_d);
                        } else if (vox_spawn_map.has(key_str)) {
                            grid->set_cell(gp, OPEN);
                            Dictionary spawn_d;
                            spawn_d["position"] = Vector3(fx, fy, fz);
                            spawn_d["palette_index"] = (int)ci;
                            spawn_d["spawn_type"] = vox_spawn_map[key_str];
                            out_spawns.push_back(spawn_d);
                        } else {
                            if (vox_inverse_density) {
                                grid->set_cell(gp, OPEN);
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
    }

    String indices_str;
    for (uint8_t idx : unique_indices) {
        if (!indices_str.is_empty()) indices_str += ", ";
        indices_str += String::num_int64(idx);
    }
    UtilityFunctions::print("UnderGenVoxStampNode: Stamped vox model '", room.vox_path, "' containing unique palette indices: [", indices_str, "]");
}

void UnderGenVoxStampNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    Dictionary context = inputs.get(0, Dictionary());
    if (context.is_empty()) return;

    Ref<DensityGrid> grid = context.get("grid", Ref<DensityGrid>());
    Array rooms_arr = context.get("rooms", Array());
    if (grid.is_null()) { outputs[0] = context; return; }

    UtilityFunctions::print("UnderGenVoxStampNode: Executing on ", rooms_arr.size(), " rooms.");

    // Re-connect and rebuild maps to ensure in-place inspector additions are captured
    for (int i = 0; i < vox_spawn_entries.size(); ++i) {
        Ref<VoxSpawnEntry> entry = vox_spawn_entries[i];
        if (entry.is_valid()) {
            if (!entry->is_connected("changed", Callable(this, "_on_entries_changed"))) {
                entry->connect("changed", Callable(this, "_on_entries_changed"));
            }
        }
    }
    for (int i = 0; i < vox_material_entries.size(); ++i) {
        Ref<VoxMaterialEntry> entry = vox_material_entries[i];
        if (entry.is_valid()) {
            if (!entry->is_connected("changed", Callable(this, "_on_entries_changed"))) {
                entry->connect("changed", Callable(this, "_on_entries_changed"));
            }
        }
    }
    _rebuild_maps_from_entries();

    _clear_vox_cache();
    std::vector<Dictionary> vox_spawns;

    for (int i = 0; i < rooms_arr.size(); ++i) {
        Dictionary r_dict = rooms_arr[i];
        String vox_path = r_dict.get("vox_path", "");
        if (vox_path.is_empty()) continue;

        ResolvedRoom room;
        room.id = r_dict.get("id", "");
        room.type = r_dict.get("type", "");
        room.position = r_dict.get("position", Vector3i());
        room.size = r_dict.get("size", Vector3i());
        room.vox_path = vox_path;
        room.exclude_from_smoothing = r_dict.get("exclude_from_smoothing", false);

        const ogt_vox_scene* scene = _load_vox(vox_path);
        if (scene) {
            Array room_connections;
            _stamp_vox(grid.ptr(), room, scene, vox_spawns, room_connections);
            r_dict["connection_points"] = room_connections;
            if (exclude_from_smoothing || room.exclude_from_smoothing) {
                r_dict["exclude_from_smoothing"] = true;
            }
        } else {
            UtilityFunctions::printerr("UnderGenVoxStampNode: Failed to load scene for vox path: '", vox_path, "'");
        }
    }

    // Pack spawns into an Array for downstream nodes
    Array spawns_array;
    for (const auto& s : vox_spawns) {
        spawns_array.append(s);
    }
    context["vox_spawns"] = spawns_array;

    UtilityFunctions::print("UnderGenVoxStampNode: Finished execution, generated ", spawns_array.size(), " spawns total.");

    outputs[0] = context;
}

} // namespace godot
