#include "undergen_astar_carver_node.h"
#include <godot_cpp/variant/utility_functions.hpp>

namespace godot {

UnderGenAStarCarverNode::UnderGenAStarCarverNode() {
    rng.instantiate();
    wobble_noise.instantiate();
}

UnderGenAStarCarverNode::~UnderGenAStarCarverNode() {}

void UnderGenAStarCarverNode::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_path_brush_min_radius", "radius"), &UnderGenAStarCarverNode::set_path_brush_min_radius);
    ClassDB::bind_method(D_METHOD("get_path_brush_min_radius"), &UnderGenAStarCarverNode::get_path_brush_min_radius);
    ClassDB::bind_method(D_METHOD("set_path_brush_max_radius", "radius"), &UnderGenAStarCarverNode::set_path_brush_max_radius);
    ClassDB::bind_method(D_METHOD("get_path_brush_max_radius"), &UnderGenAStarCarverNode::get_path_brush_max_radius);
    ClassDB::bind_method(D_METHOD("set_use_square_brush", "enabled"), &UnderGenAStarCarverNode::set_use_square_brush);
    ClassDB::bind_method(D_METHOD("get_use_square_brush"), &UnderGenAStarCarverNode::get_use_square_brush);
    ClassDB::bind_method(D_METHOD("set_vertical_movement_cost_multiplier", "mult"), &UnderGenAStarCarverNode::set_vertical_movement_cost_multiplier);
    ClassDB::bind_method(D_METHOD("get_vertical_movement_cost_multiplier"), &UnderGenAStarCarverNode::get_vertical_movement_cost_multiplier);

    ADD_PROPERTY(PropertyInfo(Variant::INT, "path_brush_min_radius"), "set_path_brush_min_radius", "get_path_brush_min_radius");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "path_brush_max_radius"), "set_path_brush_max_radius", "get_path_brush_max_radius");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "use_square_brush"), "set_use_square_brush", "get_use_square_brush");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "vertical_movement_cost_multiplier"), "set_vertical_movement_cost_multiplier", "get_vertical_movement_cost_multiplier");

    ClassDB::bind_method(D_METHOD("set_connect_from_ground_level", "enabled"), &UnderGenAStarCarverNode::set_connect_from_ground_level);
    ClassDB::bind_method(D_METHOD("get_connect_from_ground_level"), &UnderGenAStarCarverNode::get_connect_from_ground_level);
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "path_connect_from_ground_level"), "set_connect_from_ground_level", "get_connect_from_ground_level");
}

void UnderGenAStarCarverNode::set_path_brush_min_radius(int p_radius) { path_brush_min_radius = p_radius; }
int UnderGenAStarCarverNode::get_path_brush_min_radius() const { return path_brush_min_radius; }
void UnderGenAStarCarverNode::set_path_brush_max_radius(int p_radius) { path_brush_max_radius = p_radius; }
int UnderGenAStarCarverNode::get_path_brush_max_radius() const { return path_brush_max_radius; }
void UnderGenAStarCarverNode::set_use_square_brush(bool p_enabled) { use_square_brush = p_enabled; }
bool UnderGenAStarCarverNode::get_use_square_brush() const { return use_square_brush; }
void UnderGenAStarCarverNode::set_vertical_movement_cost_multiplier(float p_mult) { vertical_movement_cost_multiplier = p_mult; }
float UnderGenAStarCarverNode::get_vertical_movement_cost_multiplier() const { return vertical_movement_cost_multiplier; }
void UnderGenAStarCarverNode::set_connect_from_ground_level(bool p_enabled) { connect_from_ground_level = p_enabled; }
bool UnderGenAStarCarverNode::get_connect_from_ground_level() const { return connect_from_ground_level; }

void UnderGenAStarCarverNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    // Port 0: Generation Context (Dictionary)
    Dictionary context = inputs.get(0, Dictionary());
    if (context.is_empty()) {
        UtilityFunctions::printerr("UnderGenAStarCarverNode: Input context is empty.");
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
        UtilityFunctions::printerr("UnderGenAStarCarverNode: Grid not found in context.");
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

    // Run carving
    PathCarver carver;
    carver.path_brush_min_radius = path_brush_min_radius;
    carver.path_brush_max_radius = path_brush_max_radius;
    carver.use_square_brush = use_square_brush;
    carver.vertical_movement_cost_multiplier = vertical_movement_cost_multiplier;
    carver.connect_from_ground_level = connect_from_ground_level;
    carver.dungeon_mode = true; // Use A* corridors

    carver.create_paths_from_edges(grid.ptr(), rng.ptr(), wobble_noise, rooms, edges);

    // Write back fallback connection points to context for downstream nodes
    for (int i = 0; i < rooms_arr.size(); ++i) {
        Dictionary r_dict = rooms_arr[i];
        r_dict["connection_points"] = rooms[i].connection_points;
    }

    // Pack back to outputs
    outputs[0] = context; // Port 0: Generation Context
}

} // namespace godot
