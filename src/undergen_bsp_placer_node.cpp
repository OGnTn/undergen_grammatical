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
            if (r.vox_path.ends_with(".sdf")) {
                Ref<FileAccess> file = FileAccess::open(r.vox_path, FileAccess::READ);
                if (file.is_valid()) {
                    char magic[4];
                    file->get_buffer((uint8_t*)magic, 4);
                    if (magic[0] == 'U' && magic[1] == 'S' && magic[2] == 'D' && magic[3] == 'F') {
                        file->get_16(); // version
                        file->get_16(); // flags
                        uint32_t sx = file->get_32();
                        uint32_t sy = file->get_32();
                        uint32_t sz = file->get_32();
                        r.size = Vector3i(sx, sy, sz);
                        size_loaded = true;
                        UtilityFunctions::print("UnderGenBSPPlacerNode: Loaded size for SDF room \"", r.type, "\" from file: ", r.size);
                    } else {
                        UtilityFunctions::printerr("UnderGenBSPPlacerNode: Invalid magic header in SDF file: ", r.vox_path);
                    }
                } else {
                    UtilityFunctions::printerr("UnderGenBSPPlacerNode: Failed to open SDF file for size: ", r.vox_path);
                }
            } else {
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

    // Run Collision-Free Random Placement
    if (!processing_rooms.empty()) {
        int safe_margin_x = Math::clamp(margin_x, 1, Math::max(1, (grid_size_x - 6) / 2));
        int safe_margin_y = Math::clamp(margin_y, 1, Math::max(1, (grid_size_y - 6) / 2));
        int safe_margin_z = Math::clamp(margin_z, 1, Math::max(1, (grid_size_z - 6) / 2));

        Vector3i grid_center = grid->get_grid_dimensions() / 2;

        // Separate fixed and unfixed rooms
        std::vector<size_t> fixed_room_indices;
        std::vector<size_t> unfixed_room_indices;
        for (size_t i = 0; i < processing_rooms.size(); ++i) {
            if (processing_rooms[i].is_fixed) {
                fixed_room_indices.push_back(i);
            } else {
                unfixed_room_indices.push_back(i);
            }
        }

        // Sort unfixed rooms by volume (descending) so we place the largest/hardest-to-fit rooms first
        std::sort(unfixed_room_indices.begin(), unfixed_room_indices.end(), [&](size_t a, size_t b) {
            Vector3i sa = processing_rooms[a].size;
            Vector3i sb = processing_rooms[b].size;
            return (sa.x * sa.y * sa.z) > (sb.x * sb.y * sb.z);
        });

        // We will maintain a list of placed room indices to check for overlaps
        std::vector<size_t> placed_room_indices = fixed_room_indices;

        for (size_t idx : unfixed_room_indices) {
            ResolvedRoom &room = processing_rooms[idx];

            Vector3i min_bounds(safe_margin_x, safe_margin_y, safe_margin_z);
            Vector3i max_bounds(
                grid_size_x - safe_margin_x - room.size.x,
                grid_size_y - safe_margin_y - room.size.y,
                grid_size_z - safe_margin_z - room.size.z
            );

            // Generate candidates using dynamic spacing to prevent overlaps
            std::vector<Vector3i> candidates;
            int current_spacing = 8;

            while (candidates.empty() && current_spacing >= 0) {
                for (int attempt = 0; attempt < 1000; ++attempt) {
                    int rx = (max_bounds.x > min_bounds.x) ? rng->randi_range(min_bounds.x, max_bounds.x) : min_bounds.x;
                    int ry = (max_bounds.y > min_bounds.y) ? rng->randi_range(min_bounds.y, max_bounds.y) : min_bounds.y;
                    int rz = (max_bounds.z > min_bounds.z) ? rng->randi_range(min_bounds.z, max_bounds.z) : min_bounds.z;
                    Vector3i proposed_pos(rx, ry, rz);

                    bool overlaps = false;
                    for (size_t p_idx : placed_room_indices) {
                        const ResolvedRoom &placed = processing_rooms[p_idx];
                        bool overlap_x = (proposed_pos.x - current_spacing < placed.position.x + placed.size.x) &&
                                         (proposed_pos.x + room.size.x + current_spacing > placed.position.x);
                        bool overlap_y = (proposed_pos.y - current_spacing < placed.position.y + placed.size.y) &&
                                         (proposed_pos.y + room.size.y + current_spacing > placed.position.y);
                        bool overlap_z = (proposed_pos.z - current_spacing < placed.position.z + placed.size.z) &&
                                         (proposed_pos.z + room.size.z + current_spacing > placed.position.z);
                        if (overlap_x && overlap_y && overlap_z) {
                            overlaps = true;
                            break;
                        }
                    }

                    if (!overlaps) {
                        candidates.push_back(proposed_pos);
                        if (candidates.size() >= 30) {
                            break;
                        }
                    }
                }

                if (candidates.empty()) {
                    current_spacing -= 2;
                }
            }

            Vector3i chosen_pos;
            if (!candidates.empty()) {
                // Sort candidates by their distance to the center of the grid
                std::sort(candidates.begin(), candidates.end(), [&](const Vector3i &a, const Vector3i &b) {
                    Vector3 center_a(a.x + room.size.x / 2.0f, a.y + room.size.y / 2.0f, a.z + room.size.z / 2.0f);
                    Vector3 center_b(b.x + room.size.x / 2.0f, b.y + room.size.y / 2.0f, b.z + room.size.z / 2.0f);
                    Vector3 center_grid(grid_center.x, grid_center.y, grid_center.z);
                    return center_a.distance_to(center_grid) < center_b.distance_to(center_grid);
                });

                float t = Math::abs(spread_ratio - 0.5f) * 2.0f;
                int chosen_idx = 0;
                if (spread_ratio < 0.5f) {
                    // Blend between a random candidate and the one closest to the center (index 0)
                    int rand_idx = rng->randi() % candidates.size();
                    chosen_idx = (int)Math::round(Math::lerp((float)rand_idx, 0.0f, t));
                } else {
                    // Blend between a random candidate and the one furthest from the center (last index)
                    int rand_idx = rng->randi() % candidates.size();
                    chosen_idx = (int)Math::round(Math::lerp((float)rand_idx, (float)(candidates.size() - 1), t));
                }
                chosen_pos = candidates[chosen_idx];
            } else {
                // Complete fallback
                int rx = (max_bounds.x > min_bounds.x) ? rng->randi_range(min_bounds.x, max_bounds.x) : min_bounds.x;
                int ry = (max_bounds.y > min_bounds.y) ? rng->randi_range(min_bounds.y, max_bounds.y) : min_bounds.y;
                int rz = (max_bounds.z > min_bounds.z) ? rng->randi_range(min_bounds.z, max_bounds.z) : min_bounds.z;
                chosen_pos = Vector3i(rx, ry, rz);
                UtilityFunctions::printerr("UnderGenBSPPlacerNode: Placement fallback (overlap) for room \"", room.type, "\" (", room.id, ")");
            }

            room.position = chosen_pos;
            placed_room_indices.push_back(idx);
            UtilityFunctions::print("UnderGenBSPPlacerNode: Room \"", room.type, "\" (", room.id, ") placed at ", room.position, " size ", room.size, " (spacing used: ", current_spacing, ")");
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

Dictionary UnderGenBSPPlacerNode::get_pipeline_input_defaults(const Dictionary &global_inputs) const {
    Dictionary defaults;
    if (global_inputs.has(0)) defaults[0] = global_inputs[0];
    return defaults;
}

} // namespace godot
