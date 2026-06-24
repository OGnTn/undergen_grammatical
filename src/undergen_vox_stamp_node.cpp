#include "undergen_vox_stamp_node.h"
#include <godot_cpp/classes/file_access.hpp>
#include <godot_cpp/variant/utility_functions.hpp>
#include <godot_cpp/core/math.hpp>

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

    ADD_PROPERTY(PropertyInfo(Variant::DICTIONARY, "vox_spawn_map"), "set_vox_spawn_map", "get_vox_spawn_map");
    ADD_PROPERTY(PropertyInfo(Variant::DICTIONARY, "vox_material_map"), "set_vox_material_map", "get_vox_material_map");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "vox_inverse_density"), "set_vox_inverse_density", "get_vox_inverse_density");
}

void UnderGenVoxStampNode::set_vox_spawn_map(const Dictionary &p_map) { vox_spawn_map = p_map; }
Dictionary UnderGenVoxStampNode::get_vox_spawn_map() const { return vox_spawn_map; }
void UnderGenVoxStampNode::set_vox_material_map(const Dictionary &p_map) { vox_material_map = p_map; }
Dictionary UnderGenVoxStampNode::get_vox_material_map() const { return vox_material_map; }
void UnderGenVoxStampNode::set_vox_inverse_density(bool p_enabled) { vox_inverse_density = p_enabled; }
bool UnderGenVoxStampNode::get_vox_inverse_density() const { return vox_inverse_density; }

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
                                       std::vector<Dictionary> &out_spawns) {
    if (!scene || !grid) return;

    int gsx = grid->get_grid_size_x();
    int gsy = grid->get_grid_size_y();
    int gsz = grid->get_grid_size_z();
    Vector3i room_origin = room.position;

    constexpr float OPEN = 0.0f;
    constexpr float SOLID = 1.0f;

    for (uint32_t i = 0; i < scene->num_instances; ++i) {
        const ogt_vox_instance& inst = scene->instances[i];
        if (inst.hidden) continue;
        const ogt_vox_model* model = scene->models[inst.model_index];
        if (!model) continue;

        int offset_x = (int)inst.transform.m30;
        int offset_y = (int)inst.transform.m31;
        int offset_z = (int)inst.transform.m32;

        // Swizzle: MV(X,Y,Z) -> Godot(X, Z, Y)
        int start_grid_x = room_origin.x + offset_x;
        int start_grid_y = room_origin.y + offset_z;
        int start_grid_z = room_origin.z + offset_y;

        int lim_min_vx = Math::max(0, -start_grid_x);
        int lim_max_vx = Math::min((int)model->size_x, gsx - start_grid_x);
        int lim_min_vy = Math::max(0, -start_grid_z);
        int lim_max_vy = Math::min((int)model->size_y, gsz - start_grid_z);
        int lim_min_vz = Math::max(0, -start_grid_y);
        int lim_max_vz = Math::min((int)model->size_z, gsy - start_grid_y);

        if (lim_min_vx >= lim_max_vx || lim_min_vy >= lim_max_vy || lim_min_vz >= lim_max_vz) continue;

        for (int vx = lim_min_vx; vx < lim_max_vx; ++vx) {
            for (int vy = lim_min_vy; vy < lim_max_vy; ++vy) {
                for (int vz = lim_min_vz; vz < lim_max_vz; ++vz) {
                    uint32_t vi = vx + (vy * model->size_x) + (vz * model->size_x * model->size_y);
                    uint8_t ci = model->voxel_data[vi];

                    int fx = start_grid_x + vx;
                    int fy = start_grid_y + vz; // Model Z -> Godot Y
                    int fz = start_grid_z + vy; // Model Y -> Godot Z
                    Vector3i gp(fx, fy, fz);

                    if (ci == 0) { // Air
                        if (!vox_inverse_density) grid->set_cell(gp, OPEN);
                    } else { // Solid
                        if (vox_spawn_map.has(ci)) {
                            grid->set_cell(gp, OPEN);
                            Dictionary spawn_d;
                            spawn_d["position"] = Vector3(fx, fy, fz);
                            spawn_d["palette_index"] = (int)ci;
                            spawn_d["spawn_type"] = vox_spawn_map[ci];
                            out_spawns.push_back(spawn_d);
                        } else {
                            if (vox_inverse_density) {
                                grid->set_cell(gp, OPEN);
                            } else {
                                grid->set_cell(gp, SOLID);
                                if (vox_material_map.has(ci)) {
                                    grid->set_material_id(gp, (int)vox_material_map[ci]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void UnderGenVoxStampNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    Dictionary context = inputs.get(0, Dictionary());
    if (context.is_empty()) return;

    Ref<DensityGrid> grid = context.get("grid", Ref<DensityGrid>());
    Array rooms_arr = context.get("rooms", Array());
    if (grid.is_null()) { outputs[0] = context; return; }

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

        const ogt_vox_scene* scene = _load_vox(vox_path);
        if (scene) {
            _stamp_vox(grid.ptr(), room, scene, vox_spawns);
        }
    }

    // Pack spawns into an Array for downstream nodes
    Array spawns_array;
    for (const auto& s : vox_spawns) {
        spawns_array.append(s);
    }
    context["vox_spawns"] = spawns_array;

    outputs[0] = context;
}

} // namespace godot
