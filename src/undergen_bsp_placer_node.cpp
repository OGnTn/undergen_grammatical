#include "undergen_bsp_placer_node.h"
#include <godot_cpp/variant/utility_functions.hpp>
#include <godot_cpp/core/math.hpp>
#include <vector>
#include <map>

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

    ADD_PROPERTY(PropertyInfo(Variant::INT, "grid_size_x"), "set_grid_size_x", "get_grid_size_x");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "grid_size_y"), "set_grid_size_y", "get_grid_size_y");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "grid_size_z"), "set_grid_size_z", "get_grid_size_z");
}

void UnderGenBSPPlacerNode::set_grid_size_x(int p_x) { grid_size_x = p_x; }
int UnderGenBSPPlacerNode::get_grid_size_x() const { return grid_size_x; }
void UnderGenBSPPlacerNode::set_grid_size_y(int p_y) { grid_size_y = p_y; }
int UnderGenBSPPlacerNode::get_grid_size_y() const { return grid_size_y; }
void UnderGenBSPPlacerNode::set_grid_size_z(int p_z) { grid_size_z = p_z; }
int UnderGenBSPPlacerNode::get_grid_size_z() const { return grid_size_z; }

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
        
        Vector3i min_s = node.get("min_size", Vector3i(5, 5, 5));
        Vector3i max_s = node.get("max_size", Vector3i(10, 10, 10));
        r.size = Vector3i(
            rng->randi_range(min_s.x, max_s.x),
            rng->randi_range(min_s.y, max_s.y),
            rng->randi_range(min_s.z, max_s.z)
        );
        
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

        GridAABB root;
        root.min = Vector3i(1, 1, 1);
        root.size = Vector3i(grid_size_x - 2, grid_size_y - 2, grid_size_z - 2);

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

            // Shuffle leaves
            for (size_t i = 0; i < leaves.size(); ++i) {
                int swap_idx = rng->randi_range(0, leaves.size() - 1);
                GridAABB temp = leaves[i];
                leaves[i] = leaves[swap_idx];
                leaves[swap_idx] = temp;
            }

            for (size_t i = 0; i < processing_rooms.size(); ++i) {
                if (i >= leaves.size()) break;
                
                GridAABB leaf = leaves[i];
                ResolvedRoom &room = processing_rooms[i];
                
                if (room.is_fixed) continue;
                
                int free_x = leaf.size.x - room.size.x;
                int free_y = leaf.size.y - room.size.y;
                int free_z = leaf.size.z - room.size.z;
                
                int offset_x = (free_x > 0) ? rng->randi_range(0, free_x) : free_x / 2;
                int offset_y = (free_y > 0) ? rng->randi_range(0, free_y) : free_y / 2;
                int offset_z = (free_z > 0) ? rng->randi_range(0, free_z) : free_z / 2;
                
                room.position = leaf.min + Vector3i(offset_x, offset_y, offset_z);
            }
        }
    }

    // Stamp Generic Rooms (Pass 1 - Non-Vox rooms)
    RoomGenerator room_gen;
    int current_zone_id = 0;
    for (const auto& room : processing_rooms) {
        if (!room.vox_path.is_empty()) continue; // Skip vox rooms
        std::vector<ResolvedRoom> single_room_vec;
        single_room_vec.push_back(room);
        room_gen.create_rooms_from_data(grid.ptr(), single_room_vec, current_zone_id);
    }

    // Export Placed Rooms list back as Array of Dictionaries
    Array placed_rooms_array;
    for (const auto &room : processing_rooms) {
        Dictionary r_dict;
        r_dict["id"] = room.id;
        r_dict["type"] = room.type;
        r_dict["position"] = room.position;
        r_dict["size"] = room.size;
        r_dict["vox_path"] = room.vox_path;
        r_dict["center"] = room.center();
        placed_rooms_array.append(r_dict);
    }

    // Pack into Generation Context
    Dictionary context;
    context["grid"] = grid;
    context["rooms"] = placed_rooms_array;
    context["edges"] = edges;

    outputs[0] = context; // Port 0: Generation Context
}

} // namespace godot
