#include "undergen_bezier_carver_node.h"
#include <godot_cpp/variant/utility_functions.hpp>

namespace godot {

UnderGenBezierCarverNode::UnderGenBezierCarverNode() {
    rng.instantiate();
    rng->randomize();
    wobble_noise.instantiate();
}

UnderGenBezierCarverNode::~UnderGenBezierCarverNode() {}

void UnderGenBezierCarverNode::_bind_methods() {
    // Brush Group
    ADD_GROUP("Brush", "path_brush_");
    ClassDB::bind_method(D_METHOD("set_path_brush_min_radius", "radius"), &UnderGenBezierCarverNode::set_path_brush_min_radius);
    ClassDB::bind_method(D_METHOD("get_path_brush_min_radius"), &UnderGenBezierCarverNode::get_path_brush_min_radius);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "path_brush_min_radius"), "set_path_brush_min_radius", "get_path_brush_min_radius");

    ClassDB::bind_method(D_METHOD("set_path_brush_max_radius", "radius"), &UnderGenBezierCarverNode::set_path_brush_max_radius);
    ClassDB::bind_method(D_METHOD("get_path_brush_max_radius"), &UnderGenBezierCarverNode::get_path_brush_max_radius);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "path_brush_max_radius"), "set_path_brush_max_radius", "get_path_brush_max_radius");

    ClassDB::bind_method(D_METHOD("set_use_square_brush", "enabled"), &UnderGenBezierCarverNode::set_use_square_brush);
    ClassDB::bind_method(D_METHOD("get_use_square_brush"), &UnderGenBezierCarverNode::get_use_square_brush);
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "path_brush_square"), "set_use_square_brush", "get_use_square_brush");

    // Bezier Group
    ADD_GROUP("Bezier", "path_");
    ClassDB::bind_method(D_METHOD("set_path_segments", "segments"), &UnderGenBezierCarverNode::set_path_segments);
    ClassDB::bind_method(D_METHOD("get_path_segments"), &UnderGenBezierCarverNode::get_path_segments);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "path_segments", PROPERTY_HINT_RANGE, "1,10,1"), "set_path_segments", "get_path_segments");

    ClassDB::bind_method(D_METHOD("set_path_bend_factor", "factor"), &UnderGenBezierCarverNode::set_path_bend_factor);
    ClassDB::bind_method(D_METHOD("get_path_bend_factor"), &UnderGenBezierCarverNode::get_path_bend_factor);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "path_bend_factor", PROPERTY_HINT_RANGE, "0.0,2.0,0.05"), "set_path_bend_factor", "get_path_bend_factor");

    ClassDB::bind_method(D_METHOD("set_path_wobble_magnitude", "magnitude"), &UnderGenBezierCarverNode::set_path_wobble_magnitude);
    ClassDB::bind_method(D_METHOD("get_path_wobble_magnitude"), &UnderGenBezierCarverNode::get_path_wobble_magnitude);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "path_wobble_magnitude", PROPERTY_HINT_RANGE, "0.0,10.0,0.1"), "set_path_wobble_magnitude", "get_path_wobble_magnitude");

    ClassDB::bind_method(D_METHOD("set_path_wobble_frequency", "frequency"), &UnderGenBezierCarverNode::set_path_wobble_frequency);
    ClassDB::bind_method(D_METHOD("get_path_wobble_frequency"), &UnderGenBezierCarverNode::get_path_wobble_frequency);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "path_wobble_frequency", PROPERTY_HINT_RANGE, "0.0,1.0,0.01"), "set_path_wobble_frequency", "get_path_wobble_frequency");

    ClassDB::bind_method(D_METHOD("set_connect_from_ground_level", "enabled"), &UnderGenBezierCarverNode::set_connect_from_ground_level);
    ClassDB::bind_method(D_METHOD("get_connect_from_ground_level"), &UnderGenBezierCarverNode::get_connect_from_ground_level);
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "path_connect_from_ground_level"), "set_connect_from_ground_level", "get_connect_from_ground_level");

    // Organic Cave Shape Group
    ADD_GROUP("Cave Shape", "cave_");
    ClassDB::bind_method(D_METHOD("set_cave_ruggedness", "value"), &UnderGenBezierCarverNode::set_cave_ruggedness);
    ClassDB::bind_method(D_METHOD("get_cave_ruggedness"), &UnderGenBezierCarverNode::get_cave_ruggedness);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "cave_ruggedness", PROPERTY_HINT_RANGE, "0.0,5.0,0.05"), "set_cave_ruggedness", "get_cave_ruggedness");

    ClassDB::bind_method(D_METHOD("set_cave_floor_ruggedness", "value"), &UnderGenBezierCarverNode::set_cave_floor_ruggedness);
    ClassDB::bind_method(D_METHOD("get_cave_floor_ruggedness"), &UnderGenBezierCarverNode::get_cave_floor_ruggedness);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "cave_floor_ruggedness", PROPERTY_HINT_RANGE, "0.0,5.0,0.05"), "set_cave_floor_ruggedness", "get_cave_floor_ruggedness");

    ClassDB::bind_method(D_METHOD("set_cave_ceiling_ruggedness", "value"), &UnderGenBezierCarverNode::set_cave_ceiling_ruggedness);
    ClassDB::bind_method(D_METHOD("get_cave_ceiling_ruggedness"), &UnderGenBezierCarverNode::get_cave_ceiling_ruggedness);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "cave_ceiling_ruggedness", PROPERTY_HINT_RANGE, "0.0,5.0,0.05"), "set_cave_ceiling_ruggedness", "get_cave_ceiling_ruggedness");

    ClassDB::bind_method(D_METHOD("set_cave_width_noise", "value"), &UnderGenBezierCarverNode::set_cave_width_noise);
    ClassDB::bind_method(D_METHOD("get_cave_width_noise"), &UnderGenBezierCarverNode::get_cave_width_noise);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "cave_width_noise", PROPERTY_HINT_RANGE, "0.0,5.0,0.05"), "set_cave_width_noise", "get_cave_width_noise");

    ClassDB::bind_method(D_METHOD("set_floor_flattening", "value"), &UnderGenBezierCarverNode::set_floor_flattening);
    ClassDB::bind_method(D_METHOD("get_floor_flattening"), &UnderGenBezierCarverNode::get_floor_flattening);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "floor_flattening", PROPERTY_HINT_RANGE, "0.0,1.0,0.05"), "set_floor_flattening", "get_floor_flattening");

    ClassDB::bind_method(D_METHOD("set_overhang_openness", "value"), &UnderGenBezierCarverNode::set_overhang_openness);
    ClassDB::bind_method(D_METHOD("get_overhang_openness"), &UnderGenBezierCarverNode::get_overhang_openness);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "overhang_openness", PROPERTY_HINT_RANGE, "0.0,2.0,0.05"), "set_overhang_openness", "get_overhang_openness");
}

// Brush getters/setters
void UnderGenBezierCarverNode::set_path_brush_min_radius(int p_radius) { path_brush_min_radius = p_radius; }
int UnderGenBezierCarverNode::get_path_brush_min_radius() const { return path_brush_min_radius; }
void UnderGenBezierCarverNode::set_path_brush_max_radius(int p_radius) { path_brush_max_radius = p_radius; }
int UnderGenBezierCarverNode::get_path_brush_max_radius() const { return path_brush_max_radius; }
void UnderGenBezierCarverNode::set_use_square_brush(bool p_enabled) { use_square_brush = p_enabled; }
bool UnderGenBezierCarverNode::get_use_square_brush() const { return use_square_brush; }

// Bezier getters/setters
void UnderGenBezierCarverNode::set_path_segments(int p_segments) { path_segments = p_segments; }
int UnderGenBezierCarverNode::get_path_segments() const { return path_segments; }
void UnderGenBezierCarverNode::set_path_bend_factor(float p_factor) { path_bend_factor = p_factor; }
float UnderGenBezierCarverNode::get_path_bend_factor() const { return path_bend_factor; }
void UnderGenBezierCarverNode::set_path_wobble_magnitude(float p_magnitude) { path_wobble_magnitude = p_magnitude; }
float UnderGenBezierCarverNode::get_path_wobble_magnitude() const { return path_wobble_magnitude; }
void UnderGenBezierCarverNode::set_path_wobble_frequency(float p_frequency) { path_wobble_frequency = p_frequency; }
float UnderGenBezierCarverNode::get_path_wobble_frequency() const { return path_wobble_frequency; }
void UnderGenBezierCarverNode::set_connect_from_ground_level(bool p_enabled) { connect_from_ground_level = p_enabled; }
bool UnderGenBezierCarverNode::get_connect_from_ground_level() const { return connect_from_ground_level; }

// Organic Cave Shape getters/setters
void UnderGenBezierCarverNode::set_cave_ruggedness(float p_value) { cave_ruggedness = p_value; }
float UnderGenBezierCarverNode::get_cave_ruggedness() const { return cave_ruggedness; }
void UnderGenBezierCarverNode::set_cave_floor_ruggedness(float p_value) { cave_floor_ruggedness = p_value; }
float UnderGenBezierCarverNode::get_cave_floor_ruggedness() const { return cave_floor_ruggedness; }
void UnderGenBezierCarverNode::set_cave_ceiling_ruggedness(float p_value) { cave_ceiling_ruggedness = p_value; }
float UnderGenBezierCarverNode::get_cave_ceiling_ruggedness() const { return cave_ceiling_ruggedness; }
void UnderGenBezierCarverNode::set_cave_width_noise(float p_value) { cave_width_noise = p_value; }
float UnderGenBezierCarverNode::get_cave_width_noise() const { return cave_width_noise; }
void UnderGenBezierCarverNode::set_floor_flattening(float p_value) { floor_flattening = p_value; }
float UnderGenBezierCarverNode::get_floor_flattening() const { return floor_flattening; }
void UnderGenBezierCarverNode::set_overhang_openness(float p_value) { overhang_openness = p_value; }
float UnderGenBezierCarverNode::get_overhang_openness() const { return overhang_openness; }

void UnderGenBezierCarverNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    // Port 0: Generation Context (Dictionary)
    Dictionary context = inputs.get(0, Dictionary());
    if (context.is_empty()) {
        UtilityFunctions::printerr("UnderGenBezierCarverNode: Input context is empty.");
        return;
    }

    Ref<DensityGrid> grid = context.get("grid", Ref<DensityGrid>());
    Array rooms_arr = context.get("rooms", Array());
    Array edges_arr = context.get("edges", Array());

    if (grid.is_null()) {
        UtilityFunctions::printerr("UnderGenBezierCarverNode: Grid not found in context.");
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

    // Run bezier/organic path carving
    PathCarver carver;
    carver.path_brush_min_radius = path_brush_min_radius;
    carver.path_brush_max_radius = path_brush_max_radius;
    carver.use_square_brush = use_square_brush;
    carver.dungeon_mode = false; // Use bezier/organic paths

    // Bezier params
    carver.path_segments = path_segments;
    carver.path_bend_factor = path_bend_factor;
    carver.path_wobble_magnitude = path_wobble_magnitude;
    carver.path_wobble_frequency = path_wobble_frequency;
    carver.connect_from_ground_level = connect_from_ground_level;

    // Organic cave shape params
    carver.cave_ruggedness = cave_ruggedness;
    carver.cave_floor_ruggedness = cave_floor_ruggedness;
    carver.cave_ceiling_ruggedness = cave_ceiling_ruggedness;
    carver.cave_width_noise = cave_width_noise;
    carver.floor_flattening = floor_flattening;
    carver.overhang_openness = overhang_openness;

    carver.create_paths_from_edges(grid.ptr(), rng.ptr(), wobble_noise, rooms, edges);

    // Pack back to outputs
    outputs[0] = context; // Port 0: Generation Context
}

} // namespace godot
