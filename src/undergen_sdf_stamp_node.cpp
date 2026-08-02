#include "undergen_sdf_stamp_node.h"
#include "level_grammar_spec.h"
#include <godot_cpp/classes/file_access.hpp>
#include <godot_cpp/variant/utility_functions.hpp>
#include <godot_cpp/core/math.hpp>
#include <climits>

namespace godot {

UnderGenSdfStampNode::UnderGenSdfStampNode() {}

UnderGenSdfStampNode::~UnderGenSdfStampNode() {
    _clear_sdf_cache();
}

void UnderGenSdfStampNode::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_sdf_material_map", "map"), &UnderGenSdfStampNode::set_sdf_material_map);
    ClassDB::bind_method(D_METHOD("get_sdf_material_map"), &UnderGenSdfStampNode::get_sdf_material_map);
    ClassDB::bind_method(D_METHOD("set_blend_mode", "mode"), &UnderGenSdfStampNode::set_blend_mode);
    ClassDB::bind_method(D_METHOD("get_blend_mode"), &UnderGenSdfStampNode::get_blend_mode);
    ClassDB::bind_method(D_METHOD("set_exclude_from_smoothing", "exclude"), &UnderGenSdfStampNode::set_exclude_from_smoothing);
    ClassDB::bind_method(D_METHOD("get_exclude_from_smoothing"), &UnderGenSdfStampNode::get_exclude_from_smoothing);
    ClassDB::bind_method(D_METHOD("set_exclude_from_warping", "exclude"), &UnderGenSdfStampNode::set_exclude_from_warping);
    ClassDB::bind_method(D_METHOD("get_exclude_from_warping"), &UnderGenSdfStampNode::get_exclude_from_warping);

    ClassDB::bind_method(D_METHOD("set_sdf_material_entries", "entries"), &UnderGenSdfStampNode::set_sdf_material_entries);
    ClassDB::bind_method(D_METHOD("set_clear_room_air", "clear"), &UnderGenSdfStampNode::set_clear_room_air);
    ClassDB::bind_method(D_METHOD("get_clear_room_air"), &UnderGenSdfStampNode::get_clear_room_air);
    ClassDB::bind_method(D_METHOD("get_sdf_material_entries"), &UnderGenSdfStampNode::get_sdf_material_entries);
    ClassDB::bind_method(D_METHOD("_on_entries_changed"), &UnderGenSdfStampNode::_on_entries_changed);

    ADD_PROPERTY(PropertyInfo(Variant::DICTIONARY, "sdf_material_map", PROPERTY_HINT_NONE, "", PROPERTY_USAGE_NO_EDITOR), "set_sdf_material_map", "get_sdf_material_map");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "blend_mode", PROPERTY_HINT_ENUM, "Overwrite,MaxUnion,SubtractCarve"), "set_blend_mode", "get_blend_mode");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "exclude_from_smoothing"), "set_exclude_from_smoothing", "get_exclude_from_smoothing");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "exclude_from_warping"), "set_exclude_from_warping", "get_exclude_from_warping");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "clear_room_air"), "set_clear_room_air", "get_clear_room_air");

    ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "sdf_material_entries", PROPERTY_HINT_TYPE_STRING,
                              vformat("%d/%d:%s", Variant::OBJECT, Variant::NIL, "VoxMaterialEntry")),
                 "set_sdf_material_entries", "get_sdf_material_entries");

    BIND_ENUM_CONSTANT(MODE_OVERWRITE);
    BIND_ENUM_CONSTANT(MODE_MAX_UNION);
    BIND_ENUM_CONSTANT(MODE_SUBTRACT_CARVE);
}

void UnderGenSdfStampNode::set_sdf_material_map(const Dictionary &p_map) {
    sdf_material_map = p_map;
    
    std::map<int, Ref<VoxMaterialEntry>> existing_entries;
    for (int i = 0; i < sdf_material_entries.size(); ++i) {
        Ref<VoxMaterialEntry> entry = sdf_material_entries[i];
        if (entry.is_valid()) {
            existing_entries[entry->get_palette_index()] = entry;
        }
    }
    
    sdf_material_entries.clear();
    Array keys = p_map.keys();
    for (int i = 0; i < keys.size(); ++i) {
        Variant key = keys[i];
        int idx = 0;
        if (key.get_type() == Variant::INT) {
            idx = (int)key;
        } else if (key.get_type() == Variant::STRING) {
            idx = ((String)key).to_int();
        }
        
        Ref<VoxMaterialEntry> entry;
        auto it = existing_entries.find(idx);
        if (it != existing_entries.end()) {
            entry = it->second;
        } else {
            entry.instantiate();
            entry->set_palette_index(idx);
        }
        
        entry->set_material_id(p_map[key]);
        
        if (!entry->is_connected("changed", Callable(this, "_on_entries_changed"))) {
            entry->connect("changed", Callable(this, "_on_entries_changed"));
        }
        sdf_material_entries.append(entry);
    }
    
    emit_changed();
}

Dictionary UnderGenSdfStampNode::get_sdf_material_map() const {
    return sdf_material_map;
}

void UnderGenSdfStampNode::set_sdf_material_entries(const TypedArray<VoxMaterialEntry> &p_entries) {
    for (int i = 0; i < sdf_material_entries.size(); ++i) {
        Ref<VoxMaterialEntry> entry = sdf_material_entries[i];
        if (entry.is_valid() && entry->is_connected("changed", Callable(this, "_on_entries_changed"))) {
            entry->disconnect("changed", Callable(this, "_on_entries_changed"));
        }
    }
    
    sdf_material_entries = p_entries;
    for (int i = 0; i < sdf_material_entries.size(); ++i) {
        Ref<VoxMaterialEntry> entry = sdf_material_entries[i];
        if (entry.is_valid()) {
            if (!entry->is_connected("changed", Callable(this, "_on_entries_changed"))) {
                entry->connect("changed", Callable(this, "_on_entries_changed"));
            }
        }
    }
    
    _rebuild_maps_from_entries();
    emit_changed();
}

TypedArray<VoxMaterialEntry> UnderGenSdfStampNode::get_sdf_material_entries() const {
    return sdf_material_entries;
}

void UnderGenSdfStampNode::_rebuild_maps_from_entries() {
    sdf_material_map.clear();
    for (int i = 0; i < sdf_material_entries.size(); ++i) {
        Ref<VoxMaterialEntry> entry = sdf_material_entries[i];
        if (entry.is_valid()) {
            sdf_material_map[entry->get_palette_index()] = entry->get_material_id();
        }
    }
}

void UnderGenSdfStampNode::_on_entries_changed() {
    _rebuild_maps_from_entries();
    emit_changed();
}

void UnderGenSdfStampNode::set_blend_mode(int p_mode) { blend_mode = p_mode; emit_changed(); }
int UnderGenSdfStampNode::get_blend_mode() const { return blend_mode; }
void UnderGenSdfStampNode::set_exclude_from_smoothing(bool p_exclude) { exclude_from_smoothing = p_exclude; emit_changed(); }
bool UnderGenSdfStampNode::get_exclude_from_smoothing() const { return exclude_from_smoothing; }
void UnderGenSdfStampNode::set_exclude_from_warping(bool p_exclude) { exclude_from_warping = p_exclude; emit_changed(); }
bool UnderGenSdfStampNode::get_exclude_from_warping() const { return exclude_from_warping; }
void UnderGenSdfStampNode::set_clear_room_air(bool p_clear) { clear_room_air = p_clear; emit_changed(); }
bool UnderGenSdfStampNode::get_clear_room_air() const { return clear_room_air; }

SDFScene* UnderGenSdfStampNode::load_sdf_file(const String &path) {
    if (path.is_empty()) return nullptr;

    Ref<FileAccess> file = FileAccess::open(path, FileAccess::READ);
    if (file.is_null()) {
        UtilityFunctions::printerr("UnderGenSdfStampNode: Failed to open .sdf file: ", path);
        return nullptr;
    }

    char magic[4];
    file->get_buffer((uint8_t*)magic, 4);
    if (magic[0] != 'U' || magic[1] != 'S' || magic[2] != 'D' || magic[3] != 'F') {
        UtilityFunctions::printerr("UnderGenSdfStampNode: Invalid magic header in .sdf file: ", path);
        return nullptr;
    }

    uint16_t version = file->get_16();
    uint16_t flags = file->get_16();
    uint32_t size_x = file->get_32();
    uint32_t size_y = file->get_32();
    uint32_t size_z = file->get_32();
    float voxel_size = file->get_float();
    float pivot_x = file->get_float();
    float pivot_y = file->get_float();
    float pivot_z = file->get_float();

    SDFScene* scene = new SDFScene();
    scene->size_x = size_x;
    scene->size_y = size_y;
    scene->size_z = size_z;
    scene->voxel_size = voxel_size;
    scene->pivot = Vector3(pivot_x, pivot_y, pivot_z);

    size_t total_cells = (size_t)size_x * size_y * size_z;

    // Read Density Payload
    scene->densities.resize(total_cells);
    for (size_t i = 0; i < total_cells; ++i) {
        scene->densities[i] = file->get_float();
    }

    // Read Material Payload
    scene->materials.resize(total_cells);
    file->get_buffer(scene->materials.data(), total_cells);

    // Read Hermite Data Payload (v2)
    if (version >= 2 && (flags & 1)) {
        uint32_t num_hermite = file->get_32();
        scene->hermite_edges.resize(num_hermite);
        for (uint32_t i = 0; i < num_hermite; ++i) {
            SDFHermiteEdge &edge = scene->hermite_edges[i];
            edge.edge_type = file->get_8();
            edge.x = file->get_16();
            edge.y = file->get_16();
            edge.z = file->get_16();
            edge.t = file->get_float();
            float nx = file->get_float();
            float ny = file->get_float();
            float nz = file->get_float();
            edge.normal = Vector3(nx, ny, nz);
        }
    }

    // Read Marker Payload
    if (file->get_position() < file->get_length()) {
        uint32_t num_markers = file->get_32();
        scene->markers.resize(num_markers);
        for (uint32_t i = 0; i < num_markers; ++i) {
            SDFMarker &m = scene->markers[i];
            char name_buf[33] = {0};
            file->get_buffer((uint8_t*)name_buf, 32);
            m.name = String(name_buf);
            float px = file->get_float();
            float py = file->get_float();
            float pz = file->get_float();
            m.pos = Vector3(px, py, pz);
            float rx = file->get_float();
            float ry = file->get_float();
            float rz = file->get_float();
            float rw = file->get_float();
            m.rot = Quaternion(rx, ry, rz, rw);
        }
    }

    UtilityFunctions::print("UnderGenSdfStampNode: Loaded ", scene->markers.size(), " spawn markers from SDF '", path, "'.");

    return scene;
}

SDFScene* UnderGenSdfStampNode::_load_sdf(const String &path) {
    if (path.is_empty()) return nullptr;
    auto it = sdf_cache.find(path);
    if (it != sdf_cache.end()) return it->second;

    SDFScene* scene = load_sdf_file(path);
    if (scene) {
        sdf_cache[path] = scene;
    }
    return scene;
}

void UnderGenSdfStampNode::_clear_sdf_cache() {
    for (auto const& [path, scene] : sdf_cache) {
        delete scene;
    }
    sdf_cache.clear();
}

void UnderGenSdfStampNode::stamp_sdf_scene(DensityGrid* grid, const ResolvedRoom &room,
                                            const SDFScene* scene,
                                            const Dictionary &material_map,
                                            int blend_mode,
                                            int zone_id,
                                            std::vector<Dictionary> &out_spawns,
                                            bool clear_room_air) {
    if (!scene || !grid) return;

    UtilityFunctions::print("UnderGenSdfStampNode: Stamping SDF '", room.vox_path, "'. Grid Size: ", scene->size_x, "x", scene->size_y, "x", scene->size_z);

    Vector3i room_origin = room.position;
    int scene_size_x = (int)scene->size_x;
    int scene_size_y = (int)scene->size_y;
    int scene_size_z = (int)scene->size_z;

    // Clear room bounds to open air (0.0f) so surrounding solid bedrock doesn't clip/enclose the SDF model
    if (clear_room_air && room.size.x > 0 && room.size.y > 0 && room.size.z > 0) {
        for (int rx = 0; rx < room.size.x; ++rx) {
            for (int ry = 0; ry < room.size.y; ++ry) {
                for (int rz = 0; rz < room.size.z; ++rz) {
                    Vector3i cpos = room_origin + Vector3i(rx, ry, rz);
                    if (grid->is_valid_position(cpos)) {
                        grid->set_cell(cpos, 0.0f); // Set to open air
                        if (zone_id > 0) {
                            grid->set_zone_at(cpos, zone_id);
                        }
                    }
                }
            }
        }
    }

    float v_scale = scene->voxel_size > 0.0f ? scene->voxel_size : 1.0f;

    // Scan active solid voxels (density > 0.05f) to determine true model extent for floor alignment & centering
    int min_vx = scene_size_x, max_vx = 0;
    int min_vy = scene_size_y, max_vy = 0;
    int min_vz = scene_size_z, max_vz = 0;
    bool has_solid = false;

    for (int vx = 0; vx < scene_size_x; ++vx) {
        for (int vy = 0; vy < scene_size_y; ++vy) {
            for (int vz = 0; vz < scene_size_z; ++vz) {
                size_t idx = vx + vy * scene_size_x + vz * scene_size_x * scene_size_y;
                if (scene->densities[idx] > 0.05f) {
                    has_solid = true;
                    if (vx < min_vx) min_vx = vx;
                    if (vx > max_vx) max_vx = vx;
                    if (vy < min_vy) min_vy = vy;
                    if (vy > max_vy) max_vy = vy;
                    if (vz < min_vz) min_vz = vz;
                    if (vz > max_vz) max_vz = vz;
                }
            }
        }
    }

    if (!has_solid) {
        min_vx = 0; max_vx = scene_size_x - 1;
        min_vy = 0; max_vy = scene_size_y - 1;
        min_vz = 0; max_vz = scene_size_z - 1;
    }

    float solid_size_x = (max_vx - min_vx + 1) * v_scale;
    float solid_size_y = (max_vy - min_vy + 1) * v_scale;
    float solid_size_z = (max_vz - min_vz + 1) * v_scale;

    float offset_x = (room.size.x - solid_size_x) * 0.5f - min_vx * v_scale;
    float offset_y = 0.0f - min_vy * v_scale; // Align bottom of model to room floor
    float offset_z = (room.size.z - solid_size_z) * 0.5f - min_vz * v_scale;

    for (int vx = 0; vx < scene_size_x; ++vx) {
        for (int vy = 0; vy < scene_size_y; ++vy) {
            for (int vz = 0; vz < scene_size_z; ++vz) {
                size_t idx = vx + vy * scene_size_x + vz * scene_size_x * scene_size_y;
                float density = scene->densities[idx];
                uint8_t mat_id = scene->materials[idx];

                int fx = room_origin.x + (int)Math::round(offset_x + vx * v_scale);
                int fy = room_origin.y + (int)Math::round(offset_y + vy * v_scale);
                int fz = room_origin.z + (int)Math::round(offset_z + vz * v_scale);

                Vector3i pos(fx, fy, fz);
                if (!grid->is_valid_position(pos)) continue;

                if (blend_mode == MODE_OVERWRITE) {
                    grid->set_cell(pos, density);
                } else if (blend_mode == MODE_MAX_UNION) {
                    float current_density = grid->get_cell(pos);
                    grid->set_cell(pos, Math::max(current_density, density));
                } else if (blend_mode == MODE_SUBTRACT_CARVE) {
                    float current_density = grid->get_cell(pos);
                    grid->set_cell(pos, Math::min(current_density, 1.0f - density));
                }

                if (mat_id != 0) {
                    int mapped_mat_id = (int)mat_id;
                    Variant key = (int)mat_id;
                    if (material_map.has(key)) {
                        mapped_mat_id = (int)material_map[key];
                    }
                    grid->set_material_id(pos, mapped_mat_id);
                }

                if (zone_id > 0) {
                    grid->set_zone_at(pos, zone_id);
                }
            }
        }
    }

    // Stamp Hermite Data
    for (const SDFHermiteEdge &edge : scene->hermite_edges) {
        int fx = room_origin.x + (int)Math::round(offset_x + edge.x * v_scale);
        int fy = room_origin.y + (int)Math::round(offset_y + edge.y * v_scale);
        int fz = room_origin.z + (int)Math::round(offset_z + edge.z * v_scale);
        Vector3i edge_pos(fx, fy, fz);
        if (grid->is_valid_position(edge_pos)) {
            grid->set_hermite_edge(edge_pos, (int)edge.edge_type, edge.t, edge.normal);
        }
    }

    // Process Spawn Markers
    for (const SDFMarker &marker : scene->markers) {
        Vector3 marker_world_pos = Vector3(
            room_origin.x + offset_x + marker.pos.x * v_scale,
            room_origin.y + offset_y + marker.pos.y * v_scale,
            room_origin.z + offset_z + marker.pos.z * v_scale
        );

        Dictionary spawn;
        spawn["type"] = marker.name;
        spawn["spawn_type"] = marker.name;
        spawn["position"] = marker_world_pos;
        spawn["rotation"] = marker.rot;
        spawn["room_id"] = room.id;
        out_spawns.push_back(spawn);

        UtilityFunctions::print("UnderGenSdfStampNode: Added spawn marker '", marker.name, "' at world pos ", marker_world_pos);
    }
}

void UnderGenSdfStampNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    Dictionary context = inputs.get(0, Dictionary());
    if (context.is_empty()) return;

    Ref<DensityGrid> grid = context.get("grid", Ref<DensityGrid>());
    Array rooms_arr = context.get("rooms", Array());
    if (grid.is_null()) { outputs[0] = context; return; }

    UtilityFunctions::print("UnderGenSdfStampNode: Executing on ", rooms_arr.size(), " rooms.");

    for (int i = 0; i < sdf_material_entries.size(); ++i) {
        Ref<VoxMaterialEntry> entry = sdf_material_entries[i];
        if (entry.is_valid()) {
            if (!entry->is_connected("changed", Callable(this, "_on_entries_changed"))) {
                entry->connect("changed", Callable(this, "_on_entries_changed"));
            }
        }
    }
    _rebuild_maps_from_entries();

    _clear_sdf_cache();
    std::vector<Dictionary> sdf_spawns;

    for (int i = 0; i < rooms_arr.size(); ++i) {
        Dictionary r_dict = rooms_arr[i];
        String vox_path = r_dict.get("vox_path", "");
        if (vox_path.is_empty()) continue;
        if (!vox_path.ends_with(".sdf")) continue;

        ResolvedRoom room;
        room.id = r_dict.get("id", "");
        room.type = r_dict.get("type", "");
        room.position = r_dict.get("position", Vector3i());
        room.size = r_dict.get("size", Vector3i());
        room.vox_path = vox_path;
        room.exclude_from_smoothing = r_dict.get("exclude_from_smoothing", false);
        room.exclude_from_warping = r_dict.get("exclude_from_warping", false);

        const SDFScene* scene = _load_sdf(vox_path);
        if (scene) {
            int zone_id = grid->register_zone_name(room.type);
            stamp_sdf_scene(grid.ptr(), room, scene, sdf_material_map, blend_mode, zone_id, sdf_spawns, clear_room_air);
            if (exclude_from_smoothing || room.exclude_from_smoothing) {
                r_dict["exclude_from_smoothing"] = true;
            }
            if (exclude_from_warping || room.exclude_from_warping) {
                r_dict["exclude_from_warping"] = true;
            }
        } else {
            UtilityFunctions::printerr("UnderGenSdfStampNode: Failed to load scene for SDF path: '", vox_path, "'");
        }
    }

    Array vox_spawns_array = context.get("vox_spawns", Array());
    for (const auto &s : sdf_spawns) {
        vox_spawns_array.append(s);
    }
    context["vox_spawns"] = vox_spawns_array;

    outputs[0] = context;
}

} // namespace godot
