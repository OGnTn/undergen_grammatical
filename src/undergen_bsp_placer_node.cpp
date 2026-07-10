#include "undergen_bsp_placer_node.h"
#include <godot_cpp/variant/utility_functions.hpp>
#include <godot_cpp/core/math.hpp>
#include <godot_cpp/classes/file_access.hpp>
#include "ogt_vox.h"
#include <vector>
#include <map>
#include <algorithm>

namespace godot {

UnderGenBSPPlacerNode::UnderGenBSPPlacerNode() {
    rng.instantiate();
}

UnderGenBSPPlacerNode::~UnderGenBSPPlacerNode() {}

void UnderGenBSPPlacerNode::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_grid_size_x", "grid_size_x"), &UnderGenBSPPlacerNode::set_grid_size_x);
    ClassDB::bind_method(D_METHOD("get_grid_size_x"), &UnderGenBSPPlacerNode::get_grid_size_x);
    ClassDB::bind_method(D_METHOD("set_grid_size_y", "grid_size_y"), &UnderGenBSPPlacerNode::set_grid_size_y);
    ClassDB::bind_method(D_METHOD("get_grid_size_y"), &UnderGenBSPPlacerNode::get_grid_size_y);
    ClassDB::bind_method(D_METHOD("set_grid_size_z", "grid_size_z"), &UnderGenBSPPlacerNode::set_grid_size_z);
    ClassDB::bind_method(D_METHOD("get_grid_size_z"), &UnderGenBSPPlacerNode::get_grid_size_z);
    ClassDB::bind_method(D_METHOD("set_surface_threshold", "threshold"), &UnderGenBSPPlacerNode::set_surface_threshold);
    ClassDB::bind_method(D_METHOD("get_surface_threshold"), &UnderGenBSPPlacerNode::get_surface_threshold);
    ClassDB::bind_method(D_METHOD("set_margin_x", "margin_x"), &UnderGenBSPPlacerNode::set_margin_x);
    ClassDB::bind_method(D_METHOD("get_margin_x"), &UnderGenBSPPlacerNode::get_margin_x);
    ClassDB::bind_method(D_METHOD("set_margin_y", "margin_y"), &UnderGenBSPPlacerNode::set_margin_y);
    ClassDB::bind_method(D_METHOD("get_margin_y"), &UnderGenBSPPlacerNode::get_margin_y);
    ClassDB::bind_method(D_METHOD("set_margin_z", "margin_z"), &UnderGenBSPPlacerNode::set_margin_z);
    ClassDB::bind_method(D_METHOD("get_margin_z"), &UnderGenBSPPlacerNode::get_margin_z);
    ClassDB::bind_method(D_METHOD("set_spread_ratio", "ratio"), &UnderGenBSPPlacerNode::set_spread_ratio);
    ClassDB::bind_method(D_METHOD("get_spread_ratio"), &UnderGenBSPPlacerNode::get_spread_ratio);

    ADD_PROPERTY(PropertyInfo(Variant::INT, "grid_size_x"), "set_grid_size_x", "get_grid_size_x");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "grid_size_y"), "set_grid_size_y", "get_grid_size_y");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "grid_size_z"), "set_grid_size_z", "get_grid_size_z");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "margin_x"), "set_margin_x", "get_margin_x");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "margin_y"), "set_margin_y", "get_margin_y");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "margin_z"), "set_margin_z", "get_margin_z");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "surface_threshold"), "set_surface_threshold", "get_surface_threshold");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "spread_ratio", PROPERTY_HINT_RANGE, "0.0,1.0,0.01"), "set_spread_ratio", "get_spread_ratio");
}

void UnderGenBSPPlacerNode::set_grid_size_x(int p_x) { grid_size_x = p_x; }
int UnderGenBSPPlacerNode::get_grid_size_x() const { return grid_size_x; }
void UnderGenBSPPlacerNode::set_grid_size_y(int p_y) { grid_size_y = p_y; }
int UnderGenBSPPlacerNode::get_grid_size_y() const { return grid_size_y; }
void UnderGenBSPPlacerNode::set_grid_size_z(int p_z) { grid_size_z = p_z; }
int UnderGenBSPPlacerNode::get_grid_size_z() const { return grid_size_z; }
void UnderGenBSPPlacerNode::set_surface_threshold(float p_threshold) { surface_threshold = p_threshold; }
float UnderGenBSPPlacerNode::get_surface_threshold() const { return surface_threshold; }

void UnderGenBSPPlacerNode::set_margin_x(int p_x) { margin_x = p_x; }
int UnderGenBSPPlacerNode::get_margin_x() const { return margin_x; }
void UnderGenBSPPlacerNode::set_margin_y(int p_y) { margin_y = p_y; }
int UnderGenBSPPlacerNode::get_margin_y() const { return margin_y; }
void UnderGenBSPPlacerNode::set_margin_z(int p_z) { margin_z = p_z; }
int UnderGenBSPPlacerNode::get_margin_z() const { return margin_z; }

void UnderGenBSPPlacerNode::set_spread_ratio(float p_ratio) { spread_ratio = p_ratio; }
float UnderGenBSPPlacerNode::get_spread_ratio() const { return spread_ratio; }

void UnderGenBSPPlacerNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    // Port 0: Seed (int)
    // Port 1: Logical Graph (Dictionary)
    int64_t seed = inputs.get(0, 12345);
    Dictionary logical_graph = inputs.get(1, Dictionary());
    
    rng->set_seed(seed);

    // Initialize Grid
    Ref<DensityGrid> grid;
    grid.instantiate();
    grid->initialize_grid(grid_size_x, grid_size_y, grid_size_z, 1.0f); // Default solid
    grid->set_surface_threshold(surface_threshold);

    Array nodes = logical_graph.get("nodes", Array());
    Array edges = logical_graph.get("edges", Array());

    // Parse Nodes into internal ResolvedRoom structs
    std::vector<ResolvedRoom> processing_rooms;
    for (int i = 0; i < nodes.size(); ++i) {
        Dictionary node = nodes[i];
        ResolvedRoom r;
        r.id = node.get("id", "room");
        r.type = node.get("type", "generic");
        r.vox_path = node.get("vox_path", "");
        r.exclude_from_smoothing = node.get("exclude_from_smoothing", false);
        r.exclude_from_warping = node.get("exclude_from_warping", false);
        
        Vector3i min_s = node.get("min_size", Vector3i(5, 5, 5));
        Vector3i max_s = node.get("max_size", Vector3i(10, 10, 10));
        bool size_loaded = false;

        if (!r.vox_path.is_empty()) {
            Ref<FileAccess> file = FileAccess::open(r.vox_path, FileAccess::READ);
            if (file.is_valid()) {
                uint64_t len = file->get_length();
                PackedByteArray buf = file->get_buffer(len);
                const ogt_vox_scene* scene = ogt_vox_read_scene(buf.ptr(), (uint32_t)len);
                if (scene) {
                    int scene_min_x = INT_MAX, scene_max_x = INT_MIN;
                    int scene_min_y = INT_MAX, scene_max_y = INT_MIN;
                    int scene_min_z = INT_MAX, scene_max_z = INT_MIN;
                    bool has_visible_instances = false;

                    for (uint32_t j = 0; j < scene->num_instances; ++j) {
                        const ogt_vox_instance& inst = scene->instances[j];
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

                    if (has_visible_instances) {
                        r.size = Vector3i(
                            scene_max_x - scene_min_x,
                            scene_max_y - scene_min_y,
                            scene_max_z - scene_min_z
                        );
                        size_loaded = true;
                        UtilityFunctions::print("UnderGenBSPPlacerNode: Loaded size for vox room \"", r.type, "\" from file: ", r.size);
                    }
                    ogt_vox_destroy_scene(scene);
                } else {
                    UtilityFunctions::printerr("UnderGenBSPPlacerNode: Failed to parse vox file for size: ", r.vox_path);
                }
            } else {
                UtilityFunctions::printerr("UnderGenBSPPlacerNode: Failed to open vox file for size: ", r.vox_path);
            }
        }

        if (!size_loaded) {
            r.size = Vector3i(
                rng->randi_range(min_s.x, max_s.x),
                rng->randi_range(min_s.y, max_s.y),
                rng->randi_range(min_s.z, max_s.z)
            );
        }
        
        // Parse Constraints
        Dictionary constraints = node.get("constraints", Dictionary());
        Vector3i grid_center = grid->get_grid_dimensions() / 2;
        if (constraints.has("fixed_pos")) {
             Vector3 fix_v = constraints["fixed_pos"];
             r.is_fixed = true;
             r.fixed_position = Vector3(grid_center) + fix_v;
             r.position = Vector3i(r.fixed_position) - (r.size / 2);
        } else {
            Vector3i max_pos = grid->get_grid_dimensions() - r.size;
            r.position = Vector3i(
                rng->randi_range(0, MAX(0, max_pos.x)),
                rng->randi_range(0, MAX(0, max_pos.y)),
                rng->randi_range(0, MAX(0, max_pos.z))
            );
        }
        
        processing_rooms.push_back(r);
    }

    // Run BSP Placement
    if (!processing_rooms.empty()) {
        struct GridAABB {
            Vector3i min;
            Vector3i size;
            int volume() const { return size.x * size.y * size.z; }
        };

        // Clamp margins to leave at least 6 voxels of center space (min_child_size = 6)
        int safe_margin_x = Math::clamp(margin_x, 1, Math::max(1, (grid_size_x - 6) / 2));
        int safe_margin_y = Math::clamp(margin_y, 1, Math::max(1, (grid_size_y - 6) / 2));
        int safe_margin_z = Math::clamp(margin_z, 1, Math::max(1, (grid_size_z - 6) / 2));

        GridAABB root;
        root.min = Vector3i(safe_margin_x, safe_margin_y, safe_margin_z);
        root.size = Vector3i(grid_size_x - safe_margin_x * 2, grid_size_y - safe_margin_y * 2, grid_size_z - safe_margin_z * 2);

        if (root.size.x > 0 && root.size.y > 0 && root.size.z > 0) {
            std::vector<GridAABB> leaves;
            leaves.push_back(root);

            int rooms_needed = processing_rooms.size();
            int safety_max = 10000;
            int iterations = 0;

            while (leaves.size() < rooms_needed && iterations < safety_max) {
                iterations++;
                int best_idx = -1;
                int max_vol = -1;
                
                for (size_t i = 0; i < leaves.size(); ++i) {
                    int vol = leaves[i].volume();
                    if (vol > max_vol) {
                        max_vol = vol;
                        best_idx = i;
                    }
                }
                
                if (best_idx == -1) break;
                
                GridAABB current = leaves[best_idx];
                int axis = 0;
                if (current.size.y > current.size.x && current.size.y > current.size.z) axis = 1;
                else if (current.size.z > current.size.x) axis = 2;
                
                if (rng->randf() > 0.8) axis = rng->randi() % 3;

                int min_child_size = 6;
                int current_len = (axis == 0) ? current.size.x : (axis == 1) ? current.size.y : current.size.z;
                
                if (current_len < min_child_size * 2) {
                    break;
                }

                int split_pos = rng->randi_range(min_child_size, current_len - min_child_size);
                GridAABB c1 = current;
                GridAABB c2 = current;
                
                if (axis == 0) {
                    c1.size.x = split_pos;
                    c2.min.x += split_pos;
                    c2.size.x -= split_pos;
                } else if (axis == 1) {
                    c1.size.y = split_pos;
                    c2.min.y += split_pos;
                    c2.size.y -= split_pos;
                } else {
                    c1.size.z = split_pos;
                    c2.min.z += split_pos;
                    c2.size.z -= split_pos;
                }
                
                leaves.erase(leaves.begin() + best_idx);
                leaves.push_back(c1);
                leaves.push_back(c2);
            }

            // Separate fixed and unfixed rooms
            std::vector<size_t> unfixed_room_indices;
            for (size_t i = 0; i < processing_rooms.size(); ++i) {
                if (!processing_rooms[i].is_fixed) {
                    unfixed_room_indices.push_back(i);
                }
            }

            // Sort unfixed rooms by volume (descending)
            std::sort(unfixed_room_indices.begin(), unfixed_room_indices.end(), [&](size_t a, size_t b) {
                Vector3i sa = processing_rooms[a].size;
                Vector3i sb = processing_rooms[b].size;
                return (sa.x * sa.y * sa.z) > (sb.x * sb.y * sb.z);
            });

            // Sort leaves by volume (descending)
            std::sort(leaves.begin(), leaves.end(), [](const GridAABB &a, const GridAABB &b) {
                return a.volume() > b.volume();
            });

            Vector3i grid_center = grid->get_grid_dimensions() / 2;

            for (size_t i = 0; i < unfixed_room_indices.size(); ++i) {
                if (i >= leaves.size()) break;
                
                GridAABB leaf = leaves[i];
                ResolvedRoom &room = processing_rooms[unfixed_room_indices[i]];
                
                int free_x = leaf.size.x - room.size.x;
                int free_y = leaf.size.y - room.size.y;
                int free_z = leaf.size.z - room.size.z;
                
                int offset_x = free_x / 2;
                int offset_y = free_y / 2;
                int offset_z = free_z / 2;
                
                if (free_x > 0) {
                    int o_rand = rng->randi_range(0, free_x);
                    int o_ideal = grid_center.x - leaf.min.x - room.size.x / 2;
                    int o_close = Math::clamp(o_ideal, 0, free_x);
                    int c_min = leaf.min.x + room.size.x / 2;
                    int c_max = leaf.min.x + free_x + room.size.x / 2;
                    int o_far = (Math::abs(c_min - grid_center.x) > Math::abs(c_max - grid_center.x)) ? 0 : free_x;
                    
                    if (spread_ratio < 0.5f) {
                        float t = spread_ratio * 2.0f;
                        offset_x = (int)Math::round(Math::lerp((float)o_close, (float)o_rand, t));
                    } else {
                        float t = (spread_ratio - 0.5f) * 2.0f;
                        offset_x = (int)Math::round(Math::lerp((float)o_rand, (float)o_far, t));
                    }
                }
                
                if (free_y > 0) {
                    int o_rand = rng->randi_range(0, free_y);
                    int o_ideal = grid_center.y - leaf.min.y - room.size.y / 2;
                    int o_close = Math::clamp(o_ideal, 0, free_y);
                    int c_min = leaf.min.y + room.size.y / 2;
                    int c_max = leaf.min.y + free_y + room.size.y / 2;
                    int o_far = (Math::abs(c_min - grid_center.y) > Math::abs(c_max - grid_center.y)) ? 0 : free_y;
                    
                    if (spread_ratio < 0.5f) {
                        float t = spread_ratio * 2.0f;
                        offset_y = (int)Math::round(Math::lerp((float)o_close, (float)o_rand, t));
                    } else {
                        float t = (spread_ratio - 0.5f) * 2.0f;
                        offset_y = (int)Math::round(Math::lerp((float)o_rand, (float)o_far, t));
                    }
                }
                
                if (free_z > 0) {
                    int o_rand = rng->randi_range(0, free_z);
                    int o_ideal = grid_center.z - leaf.min.z - room.size.z / 2;
                    int o_close = Math::clamp(o_ideal, 0, free_z);
                    int c_min = leaf.min.z + room.size.z / 2;
                    int c_max = leaf.min.z + free_z + room.size.z / 2;
                    int o_far = (Math::abs(c_min - grid_center.z) > Math::abs(c_max - grid_center.z)) ? 0 : free_z;
                    
                    if (spread_ratio < 0.5f) {
                        float t = spread_ratio * 2.0f;
                        offset_z = (int)Math::round(Math::lerp((float)o_close, (float)o_rand, t));
                    } else {
                        float t = (spread_ratio - 0.5f) * 2.0f;
                        offset_z = (int)Math::round(Math::lerp((float)o_rand, (float)o_far, t));
                    }
                }
                
                room.position = leaf.min + Vector3i(offset_x, offset_y, offset_z);
            }
        }
    }

    // Stamp Generic Rooms (Pass 1 - Non-Vox rooms)
    RoomGenerator room_gen;
    int current_zone_id = 0;
    int rooms_skipped_vox = 0;
    int rooms_stamped = 0;
    for (const auto& room : processing_rooms) {
        if (!room.vox_path.is_empty()) {
            rooms_skipped_vox++;
            UtilityFunctions::print("UnderGenBSPPlacerNode: Skipping vox room \"", room.type, "\" (", room.id, ") — has vox_path: ", room.vox_path);
            continue;
        }
        std::vector<ResolvedRoom> single_room_vec;
        single_room_vec.push_back(room);
        room_gen.create_rooms_from_data(grid.ptr(), single_room_vec, current_zone_id);
        rooms_stamped++;
        UtilityFunctions::print("UnderGenBSPPlacerNode: Stamped room \"", room.type, "\" (", room.id, ") at ", room.position, " size ", room.size);
    }
    UtilityFunctions::print("UnderGenBSPPlacerNode: Room summary — ", rooms_stamped, " stamped, ", rooms_skipped_vox, " skipped (vox), ", (int)processing_rooms.size(), " total");

    // Export Placed Rooms list back as Array of Dictionaries
    Array placed_rooms_array;
    for (const auto &room : processing_rooms) {
        Dictionary r_dict;
        r_dict["id"] = room.id;
        r_dict["type"] = room.type;
        r_dict["position"] = room.position;
        r_dict["size"] = room.size;
        r_dict["vox_path"] = room.vox_path;
        r_dict["exclude_from_smoothing"] = room.exclude_from_smoothing;
        r_dict["exclude_from_warping"] = room.exclude_from_warping;
        r_dict["center"] = room.center();
        placed_rooms_array.append(r_dict);
    }

    // Pack into Generation Context
    Dictionary context;
    context["grid"] = grid;
    context["rooms"] = placed_rooms_array;
    context["edges"] = edges;
    context["seed"] = seed;

    outputs[0] = context; // Port 0: Generation Context
}

} // namespace godot
