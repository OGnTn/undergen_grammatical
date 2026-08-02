#include "undergen_modular_astar_carver_node.h"
#include <godot_cpp/classes/file_access.hpp>
#include <godot_cpp/variant/utility_functions.hpp>
#include <godot_cpp/core/math.hpp>
#include <set>
#include <algorithm>
#include <climits>

// We include ogt_vox.h without OGT_VOX_IMPLEMENTATION to avoid duplicate symbol definitions during linking
#include "ogt_vox.h"

namespace godot {

UnderGenModularAStarCarverNode::UnderGenModularAStarCarverNode() {
    rng.instantiate();
    wobble_noise.instantiate();
}

UnderGenModularAStarCarverNode::~UnderGenModularAStarCarverNode() {
    _clear_vox_cache();
}

void UnderGenModularAStarCarverNode::_bind_methods() {
    // Pathfinding Config Bindings
    ClassDB::bind_method(D_METHOD("set_path_brush_min_radius", "radius"), &UnderGenModularAStarCarverNode::set_path_brush_min_radius);
    ClassDB::bind_method(D_METHOD("get_path_brush_min_radius"), &UnderGenModularAStarCarverNode::get_path_brush_min_radius);
    ClassDB::bind_method(D_METHOD("set_path_brush_max_radius", "radius"), &UnderGenModularAStarCarverNode::set_path_brush_max_radius);
    ClassDB::bind_method(D_METHOD("get_path_brush_max_radius"), &UnderGenModularAStarCarverNode::get_path_brush_max_radius);
    ClassDB::bind_method(D_METHOD("set_use_square_brush", "enabled"), &UnderGenModularAStarCarverNode::set_use_square_brush);
    ClassDB::bind_method(D_METHOD("get_use_square_brush"), &UnderGenModularAStarCarverNode::get_use_square_brush);
    ClassDB::bind_method(D_METHOD("set_vertical_movement_cost_multiplier", "mult"), &UnderGenModularAStarCarverNode::set_vertical_movement_cost_multiplier);
    ClassDB::bind_method(D_METHOD("get_vertical_movement_cost_multiplier"), &UnderGenModularAStarCarverNode::get_vertical_movement_cost_multiplier);
    ClassDB::bind_method(D_METHOD("set_connect_from_ground_level", "enabled"), &UnderGenModularAStarCarverNode::set_connect_from_ground_level);
    ClassDB::bind_method(D_METHOD("get_connect_from_ground_level"), &UnderGenModularAStarCarverNode::get_connect_from_ground_level);

    ADD_PROPERTY(PropertyInfo(Variant::INT, "path_brush_min_radius"), "set_path_brush_min_radius", "get_path_brush_min_radius");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "path_brush_max_radius"), "set_path_brush_max_radius", "get_path_brush_max_radius");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "use_square_brush"), "set_use_square_brush", "get_use_square_brush");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "vertical_movement_cost_multiplier"), "set_vertical_movement_cost_multiplier", "get_vertical_movement_cost_multiplier");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "path_connect_from_ground_level"), "set_connect_from_ground_level", "get_connect_from_ground_level");

    // VOX Tile Paths Bindings
    ClassDB::bind_method(D_METHOD("set_doorway_vox_path", "path"), &UnderGenModularAStarCarverNode::set_doorway_vox_path);
    ClassDB::bind_method(D_METHOD("get_doorway_vox_path"), &UnderGenModularAStarCarverNode::get_doorway_vox_path);
    ClassDB::bind_method(D_METHOD("set_straight_vox_path", "path"), &UnderGenModularAStarCarverNode::set_straight_vox_path);
    ClassDB::bind_method(D_METHOD("get_straight_vox_path"), &UnderGenModularAStarCarverNode::get_straight_vox_path);
    ClassDB::bind_method(D_METHOD("set_turn_vox_path", "path"), &UnderGenModularAStarCarverNode::set_turn_vox_path);
    ClassDB::bind_method(D_METHOD("get_turn_vox_path"), &UnderGenModularAStarCarverNode::get_turn_vox_path);
    ClassDB::bind_method(D_METHOD("set_stairs_vox_path", "path"), &UnderGenModularAStarCarverNode::set_stairs_vox_path);
    ClassDB::bind_method(D_METHOD("get_stairs_vox_path"), &UnderGenModularAStarCarverNode::get_stairs_vox_path);
    ClassDB::bind_method(D_METHOD("set_t_junction_vox_path", "path"), &UnderGenModularAStarCarverNode::set_t_junction_vox_path);
    ClassDB::bind_method(D_METHOD("get_t_junction_vox_path"), &UnderGenModularAStarCarverNode::get_t_junction_vox_path);

    ADD_GROUP("Modular VOX Prefabs", "");
    ADD_PROPERTY(PropertyInfo(Variant::STRING, "doorway_vox_path", PROPERTY_HINT_FILE, "*.vox"), "set_doorway_vox_path", "get_doorway_vox_path");
    ADD_PROPERTY(PropertyInfo(Variant::STRING, "straight_vox_path", PROPERTY_HINT_FILE, "*.vox"), "set_straight_vox_path", "get_straight_vox_path");
    ADD_PROPERTY(PropertyInfo(Variant::STRING, "turn_vox_path", PROPERTY_HINT_FILE, "*.vox"), "set_turn_vox_path", "get_turn_vox_path");
    ADD_PROPERTY(PropertyInfo(Variant::STRING, "stairs_vox_path", PROPERTY_HINT_FILE, "*.vox"), "set_stairs_vox_path", "get_stairs_vox_path");
    ADD_PROPERTY(PropertyInfo(Variant::STRING, "t_junction_vox_path", PROPERTY_HINT_FILE, "*.vox"), "set_t_junction_vox_path", "get_t_junction_vox_path");

    // VOX Stamp Customization Bindings
    ClassDB::bind_method(D_METHOD("set_connection_palette_index", "palette_index"), &UnderGenModularAStarCarverNode::set_connection_palette_index);
    ClassDB::bind_method(D_METHOD("get_connection_palette_index"), &UnderGenModularAStarCarverNode::get_connection_palette_index);
    ClassDB::bind_method(D_METHOD("set_vox_inverse_density", "enabled"), &UnderGenModularAStarCarverNode::set_vox_inverse_density);
    ClassDB::bind_method(D_METHOD("get_vox_inverse_density"), &UnderGenModularAStarCarverNode::get_vox_inverse_density);
    ClassDB::bind_method(D_METHOD("set_exclude_from_smoothing", "exclude"), &UnderGenModularAStarCarverNode::set_exclude_from_smoothing);
    ClassDB::bind_method(D_METHOD("get_exclude_from_smoothing"), &UnderGenModularAStarCarverNode::get_exclude_from_smoothing);
    ClassDB::bind_method(D_METHOD("set_exclude_from_warping", "exclude"), &UnderGenModularAStarCarverNode::set_exclude_from_warping);
    ClassDB::bind_method(D_METHOD("get_exclude_from_warping"), &UnderGenModularAStarCarverNode::get_exclude_from_warping);

    ClassDB::bind_method(D_METHOD("set_vox_spawn_map", "map"), &UnderGenModularAStarCarverNode::set_vox_spawn_map);
    ClassDB::bind_method(D_METHOD("get_vox_spawn_map"), &UnderGenModularAStarCarverNode::get_vox_spawn_map);
    ClassDB::bind_method(D_METHOD("set_vox_material_map", "map"), &UnderGenModularAStarCarverNode::set_vox_material_map);
    ClassDB::bind_method(D_METHOD("get_vox_material_map"), &UnderGenModularAStarCarverNode::get_vox_material_map);
    ClassDB::bind_method(D_METHOD("set_vox_spawn_entries", "entries"), &UnderGenModularAStarCarverNode::set_vox_spawn_entries);
    ClassDB::bind_method(D_METHOD("get_vox_spawn_entries"), &UnderGenModularAStarCarverNode::get_vox_spawn_entries);
    ClassDB::bind_method(D_METHOD("set_vox_material_entries", "entries"), &UnderGenModularAStarCarverNode::set_vox_material_entries);
    ClassDB::bind_method(D_METHOD("get_vox_material_entries"), &UnderGenModularAStarCarverNode::get_vox_material_entries);

    ClassDB::bind_method(D_METHOD("_on_entries_changed"), &UnderGenModularAStarCarverNode::_on_entries_changed);

    ADD_GROUP("VOX Customizations", "");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "connection_palette_index"), "set_connection_palette_index", "get_connection_palette_index");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "vox_inverse_density"), "set_vox_inverse_density", "get_vox_inverse_density");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "exclude_from_smoothing"), "set_exclude_from_smoothing", "get_exclude_from_smoothing");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "exclude_from_warping"), "set_exclude_from_warping", "get_exclude_from_warping");

    ADD_PROPERTY(PropertyInfo(Variant::DICTIONARY, "vox_spawn_map", PROPERTY_HINT_NONE, "", PROPERTY_USAGE_NO_EDITOR), "set_vox_spawn_map", "get_vox_spawn_map");
    ADD_PROPERTY(PropertyInfo(Variant::DICTIONARY, "vox_material_map", PROPERTY_HINT_NONE, "", PROPERTY_USAGE_NO_EDITOR), "set_vox_material_map", "get_vox_material_map");

    ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "vox_spawn_entries", PROPERTY_HINT_TYPE_STRING,
                              vformat("%d/%d:%s", Variant::OBJECT, Variant::NIL, "VoxSpawnEntry")),
                 "set_vox_spawn_entries", "get_vox_spawn_entries");
    ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "vox_material_entries", PROPERTY_HINT_TYPE_STRING,
                              vformat("%d/%d:%s", Variant::OBJECT, Variant::NIL, "VoxMaterialEntry")),
                 "set_vox_material_entries", "get_vox_material_entries");
}

// Pathfinding Config Getters/Setters
void UnderGenModularAStarCarverNode::set_path_brush_min_radius(int p_radius) { path_brush_min_radius = p_radius; }
int UnderGenModularAStarCarverNode::get_path_brush_min_radius() const { return path_brush_min_radius; }
void UnderGenModularAStarCarverNode::set_path_brush_max_radius(int p_radius) { path_brush_max_radius = p_radius; }
int UnderGenModularAStarCarverNode::get_path_brush_max_radius() const { return path_brush_max_radius; }
void UnderGenModularAStarCarverNode::set_use_square_brush(bool p_enabled) { use_square_brush = p_enabled; }
bool UnderGenModularAStarCarverNode::get_use_square_brush() const { return use_square_brush; }
void UnderGenModularAStarCarverNode::set_vertical_movement_cost_multiplier(float p_mult) { vertical_movement_cost_multiplier = p_mult; }
float UnderGenModularAStarCarverNode::get_vertical_movement_cost_multiplier() const { return vertical_movement_cost_multiplier; }
void UnderGenModularAStarCarverNode::set_connect_from_ground_level(bool p_enabled) { connect_from_ground_level = p_enabled; }
bool UnderGenModularAStarCarverNode::get_connect_from_ground_level() const { return connect_from_ground_level; }

// VOX Paths Getters/Setters
void UnderGenModularAStarCarverNode::set_doorway_vox_path(const String &p_path) { doorway_vox_path = p_path; }
String UnderGenModularAStarCarverNode::get_doorway_vox_path() const { return doorway_vox_path; }
void UnderGenModularAStarCarverNode::set_straight_vox_path(const String &p_path) { straight_vox_path = p_path; }
String UnderGenModularAStarCarverNode::get_straight_vox_path() const { return straight_vox_path; }
void UnderGenModularAStarCarverNode::set_turn_vox_path(const String &p_path) { turn_vox_path = p_path; }
String UnderGenModularAStarCarverNode::get_turn_vox_path() const { return turn_vox_path; }
void UnderGenModularAStarCarverNode::set_stairs_vox_path(const String &p_path) { stairs_vox_path = p_path; }
String UnderGenModularAStarCarverNode::get_stairs_vox_path() const { return stairs_vox_path; }
void UnderGenModularAStarCarverNode::set_t_junction_vox_path(const String &p_path) { t_junction_vox_path = p_path; }
String UnderGenModularAStarCarverNode::get_t_junction_vox_path() const { return t_junction_vox_path; }

// VOX Stamp Customization Getters/Setters
void UnderGenModularAStarCarverNode::set_connection_palette_index(int p_index) { connection_palette_index = p_index; }
int UnderGenModularAStarCarverNode::get_connection_palette_index() const { return connection_palette_index; }
void UnderGenModularAStarCarverNode::set_vox_inverse_density(bool p_enabled) { vox_inverse_density = p_enabled; }
bool UnderGenModularAStarCarverNode::get_vox_inverse_density() const { return vox_inverse_density; }
void UnderGenModularAStarCarverNode::set_exclude_from_smoothing(bool p_exclude) { exclude_from_smoothing = p_exclude; }
bool UnderGenModularAStarCarverNode::get_exclude_from_smoothing() const { return exclude_from_smoothing; }
void UnderGenModularAStarCarverNode::set_exclude_from_warping(bool p_exclude) { exclude_from_warping = p_exclude; }
bool UnderGenModularAStarCarverNode::get_exclude_from_warping() const { return exclude_from_warping; }

void UnderGenModularAStarCarverNode::set_vox_spawn_map(const Dictionary &p_map) {
    vox_spawn_map = p_map;
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
        int idx = (key.get_type() == Variant::INT) ? (int)key : ((String)key).to_int();
        entry->set_palette_index(idx);
        entry->set_spawn_type(p_map[key]);
        entry->connect("changed", Callable(this, "_on_entries_changed"));
        vox_spawn_entries.append(entry);
    }
}
Dictionary UnderGenModularAStarCarverNode::get_vox_spawn_map() const { return vox_spawn_map; }

void UnderGenModularAStarCarverNode::set_vox_material_map(const Dictionary &p_map) {
    vox_material_map = p_map;
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
        int idx = (key.get_type() == Variant::INT) ? (int)key : ((String)key).to_int();
        entry->set_palette_index(idx);
        entry->set_material_id(p_map[key]);
        entry->connect("changed", Callable(this, "_on_entries_changed"));
        vox_material_entries.append(entry);
    }
}
Dictionary UnderGenModularAStarCarverNode::get_vox_material_map() const { return vox_material_map; }

void UnderGenModularAStarCarverNode::set_vox_spawn_entries(const TypedArray<VoxSpawnEntry> &p_entries) {
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
}
TypedArray<VoxSpawnEntry> UnderGenModularAStarCarverNode::get_vox_spawn_entries() const { return vox_spawn_entries; }

void UnderGenModularAStarCarverNode::set_vox_material_entries(const TypedArray<VoxMaterialEntry> &p_entries) {
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
}
TypedArray<VoxMaterialEntry> UnderGenModularAStarCarverNode::get_vox_material_entries() const { return vox_material_entries; }

void UnderGenModularAStarCarverNode::_rebuild_maps_from_entries() {
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

void UnderGenModularAStarCarverNode::_on_entries_changed() {
    _rebuild_maps_from_entries();
}

const ogt_vox_scene* UnderGenModularAStarCarverNode::_load_vox(const String &path) {
    if (path.is_empty()) return nullptr;
    auto it = vox_cache.find(path);
    if (it != vox_cache.end()) return it->second;

    Ref<FileAccess> file = FileAccess::open(path, FileAccess::READ);
    if (file.is_null()) {
        UtilityFunctions::printerr("UnderGenModularAStarCarverNode: Failed to open vox file: ", path);
        return nullptr;
    }
    uint64_t len = file->get_length();
    PackedByteArray buf = file->get_buffer(len);
    const ogt_vox_scene* scene = ogt_vox_read_scene(buf.ptr(), (uint32_t)len);
    if (!scene) {
        UtilityFunctions::printerr("UnderGenModularAStarCarverNode: Failed to parse vox file: ", path);
        return nullptr;
    }
    vox_cache[path] = scene;
    return scene;
}

void UnderGenModularAStarCarverNode::_clear_vox_cache() {
    for (auto const& [path, scene] : vox_cache) {
        ogt_vox_destroy_scene(scene);
    }
    vox_cache.clear();
}

void UnderGenModularAStarCarverNode::_stamp_vox_tile(DensityGrid* grid, const Vector3i &grid_pos, const ogt_vox_scene* scene, int zone_id, int rotation_idx, std::vector<Dictionary> &out_spawns, Array &out_rooms) {
    if (!scene || !grid) return;

    int gsx = grid->get_grid_size_x();
    int gsy = grid->get_grid_size_y();
    int gsz = grid->get_grid_size_z();

    constexpr float OPEN = 0.0f;
    constexpr float SOLID = 1.0f;

    // Scan instances to get bounding box in Godot space
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
        int inst_min_y = offset_z; // MV Z -> Godot Y
        int inst_max_y = offset_z + model->size_z;
        int inst_min_z = offset_y; // MV Y -> Godot Z
        int inst_max_z = offset_y + model->size_y;

        if (inst_min_x < scene_min_x) scene_min_x = inst_min_x;
        if (inst_max_x > scene_max_x) scene_max_x = inst_max_x;
        if (inst_min_y < scene_min_y) scene_min_y = inst_min_y;
        if (inst_max_y > scene_max_y) scene_max_y = inst_max_y;
        if (inst_min_z < scene_min_z) scene_min_z = inst_min_z;
        if (inst_max_z > scene_max_z) scene_max_z = inst_max_z;
    }

    if (!has_visible_instances) return;

    int scene_size_x = scene_max_x - scene_min_x;
    int scene_size_y = scene_max_y - scene_min_y;
    int scene_size_z = scene_max_z - scene_min_z;

    if (exclude_from_smoothing || exclude_from_warping) {
        int sz_x = (rotation_idx % 2 == 0) ? scene_size_x : scene_size_z;
        int sz_z = (rotation_idx % 2 == 0) ? scene_size_z : scene_size_x;

        Dictionary r_dict;
        r_dict["id"] = "corr_stamp_" + String::num_int64(grid_pos.x) + "_" + String::num_int64(grid_pos.y) + "_" + String::num_int64(grid_pos.z);
        r_dict["type"] = "corridor";
        r_dict["position"] = Vector3i(grid_pos.x - sz_x / 2, grid_pos.y, grid_pos.z - sz_z / 2);
        r_dict["size"] = Vector3i(sz_x, scene_size_y, sz_z);
        r_dict["exclude_from_smoothing"] = exclude_from_smoothing;
        r_dict["exclude_from_warping"] = exclude_from_warping;
        out_rooms.append(r_dict);
    }

    // Stamp each voxel with integer-based Y-rotation applied relative to center
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

                    int rx = offset_x + vx;
                    int ry = offset_z + vz;
                    int rz = offset_y + vy;

                    int lx = rx - scene_min_x;
                    int ly = ry - scene_min_y;
                    int lz = rz - scene_min_z;

                    int ox = 0;
                    int oy = ly; // Align bottom of tile with grid_pos.y
                    int oz = 0;

                    switch (rotation_idx) {
                        case 0:
                            ox = lx - scene_size_x / 2;
                            oz = lz - scene_size_z / 2;
                            break;
                        case 1: // 90 CW
                            ox = (scene_size_z - 1 - lz) - scene_size_z / 2;
                            oz = lx - scene_size_x / 2;
                            break;
                        case 2: // 180
                            ox = (scene_size_x - 1 - lx) - scene_size_x / 2;
                            oz = (scene_size_z - 1 - lz) - scene_size_z / 2;
                            break;
                        case 3: // 270 CW
                            ox = lz - scene_size_z / 2;
                            oz = (scene_size_x - 1 - lx) - scene_size_x / 2;
                            break;
                    }

                    int fx = grid_pos.x + ox;
                    int fy = grid_pos.y + oy;
                    int fz = grid_pos.z + oz;

                    if (fx < 0 || fx >= gsx || fy < 0 || fy >= gsy || fz < 0 || fz >= gsz) {
                        continue;
                    }

                    Vector3i gp(fx, fy, fz);
                    if (zone_id > 0) {
                        grid->set_zone_at(gp, zone_id);
                    }

                    if (ci == 0) {
                        if (!vox_inverse_density) grid->set_cell(gp, OPEN);
                    } else {
                        Variant key_int = (int)ci;
                        Variant key_str = String::num_int64(ci);
                        if (ci == connection_palette_index) {
                            grid->set_cell(gp, OPEN);
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
}

void UnderGenModularAStarCarverNode::_stamp_vox_tile_fallback(DensityGrid* grid, const Vector3i& center, int zone_id) {
    int r_low = path_brush_min_radius;
    int r_high = path_brush_max_radius;
    if (r_high < r_low) r_high = r_low;
    if (r_low < 0) r_low = 0;

    int radius = (r_low == r_high) ? r_low : rng->randi_range(r_low, r_high);
    if (radius <= 0) return;

    int gsx = grid->get_grid_size_x();
    int gsy = grid->get_grid_size_y();
    int gsz = grid->get_grid_size_z();
    constexpr float OPEN = 0.0f;

    for (int i = -radius; i <= radius; ++i) {
        for (int j = 0; j <= radius * 2; ++j) { 
            for (int k = -radius; k <= radius; ++k) {
                Vector3i current_pos = center + Vector3i(i, j, k);
                
                if (current_pos.x > 0 && current_pos.x < gsx - 1 &&
                    current_pos.y > 0 && current_pos.y < gsy - 1 &&
                    current_pos.z > 0 && current_pos.z < gsz - 1) {
                    
                    grid->set_cell(current_pos, OPEN);
                    if (zone_id > 0 && grid->is_valid_position(current_pos)) {
                        if (grid->get_zone_at(current_pos) == 0) {
                            grid->set_zone_at(current_pos, zone_id);
                        }
                    }
                }
            }
        }
    }
}

void UnderGenModularAStarCarverNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    Dictionary context = inputs.get(0, Dictionary());
    if (context.is_empty()) {
        UtilityFunctions::printerr("UnderGenModularAStarCarverNode: Input context is empty.");
        return;
    }

    int64_t seed = context.get("seed", (int64_t)12345);
    rng->set_seed(seed);
    if (wobble_noise.is_valid()) {
        wobble_noise->set_seed((int)seed);
    }

    Ref<DensityGrid> grid = context.get("grid", Ref<DensityGrid>());
    Array rooms_arr = context.get("rooms", Array());
    Array edges_arr = context.get("edges", Array());

    if (grid.is_null()) {
        UtilityFunctions::printerr("UnderGenModularAStarCarverNode: Grid not found in context.");
        return;
    }

    // Reconstruct ResolvedRoom vector
    std::vector<ResolvedRoom> rooms;
    std::map<String, int> id_to_index;
    for (int i = 0; i < rooms_arr.size(); ++i) {
        Dictionary r_dict = rooms_arr[i];
        ResolvedRoom r;
        r.id = r_dict.get("id", "");
        r.type = r_dict.get("type", "");
        r.position = r_dict.get("position", Vector3i());
        r.size = r_dict.get("size", Vector3i());
        r.vox_path = r_dict.get("vox_path", "");
        r.connection_points = r_dict.get("connection_points", Array());
        rooms.push_back(r);
        id_to_index[r.id] = i;
    }

    // Reconstruct ResolvedEdge vector
    std::vector<ResolvedEdge> edges;
    for (int i = 0; i < edges_arr.size(); ++i) {
        Dictionary e_dict = edges_arr[i];
        String from = e_dict.get("from", "");
        String to = e_dict.get("to", "");
        if (id_to_index.count(from) && id_to_index.count(to)) {
            ResolvedEdge e;
            e.from_index = id_to_index[from];
            e.to_index = id_to_index[to];
            e.type = e_dict.get("type", "corridor");
            edges.push_back(e);
        }
    }

    // Resolve paths via PathCarver
    PathCarver carver;
    carver.path_brush_min_radius = path_brush_min_radius;
    carver.path_brush_max_radius = path_brush_max_radius;
    carver.use_square_brush = use_square_brush;
    carver.vertical_movement_cost_multiplier = vertical_movement_cost_multiplier;
    carver.connect_from_ground_level = connect_from_ground_level;
    carver.dungeon_mode = true;

    std::vector<std::vector<Vector3i>> paths = carver.get_dungeon_paths(grid.ptr(), rng.ptr(), rooms, edges);

    // Save fallback connection points back to context
    for (int i = 0; i < rooms_arr.size(); ++i) {
        Dictionary r_dict = rooms_arr[i];
        r_dict["connection_points"] = rooms[i].connection_points;
    }

    // Build set of doorway coordinates
    std::set<Vector3i> doorway_points;
    for (const auto& room : rooms) {
        for (int j = 0; j < room.connection_points.size(); ++j) {
            doorway_points.insert(Vector3i(room.connection_points[j]));
        }
    }

    // Build the combined corridor graph and record edge types
    std::map<Vector3i, std::set<Vector3i>> corridor_graph;
    std::map<Vector3i, String> coord_to_zone_type;
    
    for (size_t path_idx = 0; path_idx < paths.size(); ++path_idx) {
        const auto& path = paths[path_idx];
        String edge_type = (path_idx < edges.size()) ? edges[path_idx].type : "corridor";
        
        for (size_t i = 0; i < path.size(); ++i) {
            Vector3i curr = path[i];
            coord_to_zone_type[curr] = edge_type;
            
            if (i > 0) {
                Vector3i prev = path[i - 1];
                corridor_graph[curr].insert(prev);
                corridor_graph[prev].insert(curr);
            }
            if (i < path.size() - 1) {
                Vector3i next = path[i + 1];
                corridor_graph[curr].insert(next);
                corridor_graph[next].insert(curr);
            }
        }
    }

    // Load VOX models
    _clear_vox_cache();
    _rebuild_maps_from_entries();

    const ogt_vox_scene* doorway_scene = _load_vox(doorway_vox_path);
    const ogt_vox_scene* straight_scene = _load_vox(straight_vox_path);
    const ogt_vox_scene* turn_scene = _load_vox(turn_vox_path);
    const ogt_vox_scene* stairs_scene = _load_vox(stairs_vox_path);
    const ogt_vox_scene* t_junction_scene = _load_vox(t_junction_vox_path);

    std::vector<Dictionary> out_spawns;

    auto get_rotation_index = [](const Vector3i& d) -> int {
        if (d.x > 0) return 0; // East (+X)
        if (d.z > 0) return 1; // South (+Z)
        if (d.x < 0) return 2; // West (-X)
        if (d.z < 0) return 3; // North (-Z)
        return 0;
    };

    // Classify and Stamp each node
    for (const auto& [pos, neighbors] : corridor_graph) {
        String zone_name = coord_to_zone_type.count(pos) ? coord_to_zone_type[pos] : "corridor";
        int zone_id = grid->register_zone_name(zone_name);

        // 1. Check if Doorway
        if (doorway_points.count(pos)) {
            if (doorway_scene) {
                Vector3i d_out(0, 0, 0);
                if (neighbors.size() > 0) {
                    d_out = *neighbors.begin() - pos;
                }
                int rot = get_rotation_index(d_out);
                _stamp_vox_tile(grid.ptr(), pos, doorway_scene, zone_id, rot, out_spawns, rooms_arr);
            } else {
                _stamp_vox_tile_fallback(grid.ptr(), pos, zone_id);
            }
            continue;
        }

        // 2. Classify based on neighbors
        if (neighbors.size() == 1) {
            // End of path corridor tile
            if (straight_scene) {
                Vector3i d = *neighbors.begin() - pos;
                int rot = get_rotation_index(d);
                _stamp_vox_tile(grid.ptr(), pos, straight_scene, zone_id, rot, out_spawns, rooms_arr);
            } else {
                _stamp_vox_tile_fallback(grid.ptr(), pos, zone_id);
            }
        }
        else if (neighbors.size() == 2) {
            auto it = neighbors.begin();
            Vector3i n1 = *it;
            Vector3i n2 = *(++it);

            Vector3i d1 = n1 - pos;
            Vector3i d2 = n2 - pos;

            // Stair segment
            if (d1.y != 0 || d2.y != 0) {
                if (stairs_scene) {
                    Vector3i d_flat = (d1.y != 0) ? Vector3i(d2.x, 0, d2.z) : Vector3i(d1.x, 0, d1.z);
                    int rot = get_rotation_index(d_flat);
                    int y_change = (d1.y != 0) ? -d1.y : d2.y;
                    if (y_change < 0) {
                        rot = (rot + 2) % 4; // Face opposite direction for downhill stair
                    }
                    _stamp_vox_tile(grid.ptr(), pos, stairs_scene, zone_id, rot, out_spawns, rooms_arr);
                } else {
                    _stamp_vox_tile_fallback(grid.ptr(), pos, zone_id);
                }
            }
            // Straight segment
            else if (d1 + d2 == Vector3i(0, 0, 0)) {
                if (straight_scene) {
                    int rot = get_rotation_index(d1);
                    _stamp_vox_tile(grid.ptr(), pos, straight_scene, zone_id, rot, out_spawns, rooms_arr);
                } else {
                    _stamp_vox_tile_fallback(grid.ptr(), pos, zone_id);
                }
            }
            // Turn segment
            else {
                if (turn_scene) {
                    int rot = 0;
                    if ((d1.z < 0 || d2.z < 0) && (d1.x > 0 || d2.x > 0)) rot = 0;      // North & East
                    else if ((d1.x > 0 || d2.x > 0) && (d1.z > 0 || d2.z > 0)) rot = 1; // East & South
                    else if ((d1.z > 0 || d2.z > 0) && (d1.x < 0 || d2.x < 0)) rot = 2; // South & West
                    else if ((d1.x < 0 || d2.x < 0) && (d1.z < 0 || d2.z < 0)) rot = 3; // West & North
                    _stamp_vox_tile(grid.ptr(), pos, turn_scene, zone_id, rot, out_spawns, rooms_arr);
                } else {
                    _stamp_vox_tile_fallback(grid.ptr(), pos, zone_id);
                }
            }
        }
        else if (neighbors.size() == 3) {
            // T-Junction
            if (t_junction_scene) {
                auto it = neighbors.begin();
                Vector3i n1 = *it;
                Vector3i n2 = *(++it);
                Vector3i n3 = *(++it);

                Vector3i d1 = n1 - pos;
                Vector3i d2 = n2 - pos;
                Vector3i d3 = n3 - pos;

                Vector3i d_leg(0, 0, 0);
                if (d1 + d2 == Vector3i(0, 0, 0)) d_leg = d3;
                else if (d1 + d3 == Vector3i(0, 0, 0)) d_leg = d2;
                else if (d2 + d3 == Vector3i(0, 0, 0)) d_leg = d1;

                int rot = get_rotation_index(d_leg);
                _stamp_vox_tile(grid.ptr(), pos, t_junction_scene, zone_id, rot, out_spawns, rooms_arr);
            } else {
                _stamp_vox_tile_fallback(grid.ptr(), pos, zone_id);
            }
        }
        else {
            // 4+ neighbors (Cross junction or overlap)
            _stamp_vox_tile_fallback(grid.ptr(), pos, zone_id);
        }
    }

    _clear_vox_cache();

    // Pack spawned items to context
    Array context_spawns = context.get("vox_spawns", Array());
    for (const auto& s : out_spawns) {
        context_spawns.append(s);
    }
    context["vox_spawns"] = context_spawns;
    context["rooms"] = rooms_arr;

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

    outputs[0] = context;
}

} // namespace godot
