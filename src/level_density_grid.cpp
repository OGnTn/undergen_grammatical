// level_density_grid.cpp
#include "level_density_grid.h"

#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/variant/utility_functions.hpp>
#include <godot_cpp/core/math.hpp> // For clamp, min, max, etc.
#include <godot_cpp/classes/global_constants.hpp> // For Error enum if needed
#include <godot_cpp/classes/node.hpp> // For _Ready() if needed, although Resources don't have it
#include <map> // For std::map if needed for A* point map instead of Dictionary
#include <vector> // Added for temporary grid in low_pass and path nodes
#include <queue>
#include <unordered_map>
#include <cmath>
#include <unordered_set>
#include <godot_cpp/classes/time.hpp>

namespace godot {

LevelDensityGrid::LevelDensityGrid() {
    // Initialize Ref<> types
    noise_generator.instantiate(); // Create a default instance
    rng.instantiate();
    rng->randomize(); // Seed the random number generator

    // Set default noise properties on the created instance
    if (noise_generator.is_valid()) {
        noise_generator->set_noise_type(FastNoiseLite::TYPE_PERLIN);
        noise_generator->set_seed(noise_seed);
        noise_generator->set_frequency(noise_frequency);
    }
}

LevelDensityGrid::~LevelDensityGrid() {
    // Destructor logic if needed (Refs are automatically handled)
}

void LevelDensityGrid::_bind_methods() {
    // --- Main Generation Function ---
    ClassDB::bind_method(D_METHOD("generate_level_data", "world_grid_dimensions", "voxel_size", "seed"), &LevelDensityGrid::generate_level_data);
    ClassDB::bind_method(D_METHOD("generate_from_graph", "node_data", "edge_data", "voxel_size", "seed"), &LevelDensityGrid::generate_from_graph);
    ClassDB::bind_method(D_METHOD("_mark_brush", "center", "radius_lower", "radius_higher", "value"), &LevelDensityGrid::_mark_brush);
    // --- Getters for Runtime Data ---
    ClassDB::bind_method(D_METHOD("get_calculated_spawn_position"), &LevelDensityGrid::get_calculated_spawn_position);
    ClassDB::bind_method(D_METHOD("get_calculated_end_position"), &LevelDensityGrid::get_calculated_end_position);
    ClassDB::bind_method(D_METHOD("get_surface_normals"), &LevelDensityGrid::get_surface_normals);
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3, "calculated_spawn_position", PROPERTY_HINT_NONE, "", PROPERTY_USAGE_READ_ONLY), "", "get_calculated_spawn_position");
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3, "calculated_end_position", PROPERTY_HINT_NONE, "", PROPERTY_USAGE_READ_ONLY), "", "get_calculated_end_position");
    ADD_PROPERTY(PropertyInfo(Variant::DICTIONARY, "surface_normals", PROPERTY_HINT_NONE, "", PROPERTY_USAGE_READ_ONLY), "", "get_surface_normals");


    // --- Property Getters/Setters ---
    // Noise Group
    ADD_GROUP("Noise", "noise_");
    ClassDB::bind_method(D_METHOD("set_noise_scale", "scale"), &LevelDensityGrid::set_noise_scale);
    ClassDB::bind_method(D_METHOD("get_noise_scale"), &LevelDensityGrid::get_noise_scale);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "noise_scale"), "set_noise_scale", "get_noise_scale");

    ClassDB::bind_method(D_METHOD("set_noise_intensity", "intensity"), &LevelDensityGrid::set_noise_intensity);
    ClassDB::bind_method(D_METHOD("get_noise_intensity"), &LevelDensityGrid::get_noise_intensity);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "noise_intensity"), "set_noise_intensity", "get_noise_intensity");

    ClassDB::bind_method(D_METHOD("set_noise_seed", "seed"), &LevelDensityGrid::set_noise_seed);
    ClassDB::bind_method(D_METHOD("get_noise_seed"), &LevelDensityGrid::get_noise_seed);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "noise_seed"), "set_noise_seed", "get_noise_seed");

    ClassDB::bind_method(D_METHOD("set_noise_frequency", "freq"), &LevelDensityGrid::set_noise_frequency);
    ClassDB::bind_method(D_METHOD("get_noise_frequency"), &LevelDensityGrid::get_noise_frequency);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "noise_frequency"), "set_noise_frequency", "get_noise_frequency");

    ClassDB::bind_method(D_METHOD("set_noise_generator", "noise"), &LevelDensityGrid::set_noise_generator);
    ClassDB::bind_method(D_METHOD("get_noise_generator"), &LevelDensityGrid::get_noise_generator);
    ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "noise_generator", PROPERTY_HINT_RESOURCE_TYPE, "FastNoiseLite"), "set_noise_generator", "get_noise_generator");


    // Rooms Group
    ADD_GROUP("Rooms", "room_");
    ClassDB::bind_method(D_METHOD("set_room_count", "count"), &LevelDensityGrid::set_room_count);
    ClassDB::bind_method(D_METHOD("get_room_count"), &LevelDensityGrid::get_room_count);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "room_count"), "set_room_count", "get_room_count");

    ClassDB::bind_method(D_METHOD("set_min_room_size", "size"), &LevelDensityGrid::set_min_room_size);
    ClassDB::bind_method(D_METHOD("get_min_room_size"), &LevelDensityGrid::get_min_room_size);
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3I, "room_min_size"), "set_min_room_size", "get_min_room_size");

    ClassDB::bind_method(D_METHOD("set_max_room_size", "size"), &LevelDensityGrid::set_max_room_size);
    ClassDB::bind_method(D_METHOD("get_max_room_size"), &LevelDensityGrid::get_max_room_size);
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3I, "room_max_size"), "set_max_room_size", "get_max_room_size");

    ClassDB::bind_method(D_METHOD("set_max_placement_tries", "tries"), &LevelDensityGrid::set_max_placement_tries);
    ClassDB::bind_method(D_METHOD("get_max_placement_tries"), &LevelDensityGrid::get_max_placement_tries);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "room_max_placement_tries"), "set_max_placement_tries", "get_max_placement_tries");


    // Paths Group
    ADD_GROUP("Paths", "path_");
    ClassDB::bind_method(D_METHOD("set_dungeon_mode", "enabled"), &LevelDensityGrid::set_dungeon_mode);
    ClassDB::bind_method(D_METHOD("get_dungeon_mode"), &LevelDensityGrid::get_dungeon_mode);
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "path_dungeon_mode"), "set_dungeon_mode", "get_dungeon_mode");

    ClassDB::bind_method(D_METHOD("set_path_brush_min_radius", "radius"), &LevelDensityGrid::set_path_brush_min_radius);
    ClassDB::bind_method(D_METHOD("get_path_brush_min_radius"), &LevelDensityGrid::get_path_brush_min_radius);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "path_brush_min_radius"), "set_path_brush_min_radius", "get_path_brush_min_radius");

    ClassDB::bind_method(D_METHOD("set_path_brush_max_radius", "radius"), &LevelDensityGrid::set_path_brush_max_radius);
    ClassDB::bind_method(D_METHOD("get_path_brush_max_radius"), &LevelDensityGrid::get_path_brush_max_radius);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "path_brush_max_radius"), "set_path_brush_max_radius", "get_path_brush_max_radius");

    ClassDB::bind_method(D_METHOD("set_use_square_brush", "enabled"), &LevelDensityGrid::set_use_square_brush);
    ClassDB::bind_method(D_METHOD("get_use_square_brush"), &LevelDensityGrid::get_use_square_brush);
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "path_use_square_brush"), "set_use_square_brush", "get_use_square_brush");

    ClassDB::bind_method(D_METHOD("set_vertical_movement_cost_multiplier", "mult"), &LevelDensityGrid::set_vertical_movement_cost_multiplier);
    ClassDB::bind_method(D_METHOD("get_vertical_movement_cost_multiplier"), &LevelDensityGrid::get_vertical_movement_cost_multiplier);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "path_vertical_movement_cost_multiplier"), "set_vertical_movement_cost_multiplier", "get_vertical_movement_cost_multiplier");
    
    // Bezier Path Bindings
    ClassDB::bind_method(D_METHOD("set_path_segments", "segments"), &LevelDensityGrid::set_path_segments);
    ClassDB::bind_method(D_METHOD("get_path_segments"), &LevelDensityGrid::get_path_segments);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "path_segments", PROPERTY_HINT_RANGE, "1,10,1"), "set_path_segments", "get_path_segments");

    ClassDB::bind_method(D_METHOD("set_path_bend_factor", "factor"), &LevelDensityGrid::set_path_bend_factor);
    ClassDB::bind_method(D_METHOD("get_path_bend_factor"), &LevelDensityGrid::get_path_bend_factor);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "path_bend_factor", PROPERTY_HINT_RANGE, "0.0,2.0,0.05"), "set_path_bend_factor", "get_path_bend_factor");
    
    ClassDB::bind_method(D_METHOD("set_path_wobble_magnitude", "magnitude"), &LevelDensityGrid::set_path_wobble_magnitude);
    ClassDB::bind_method(D_METHOD("get_path_wobble_magnitude"), &LevelDensityGrid::get_path_wobble_magnitude);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "path_wobble_magnitude", PROPERTY_HINT_RANGE, "0.0,10.0,0.1"), "set_path_wobble_magnitude", "get_path_wobble_magnitude");
    
    ClassDB::bind_method(D_METHOD("set_path_wobble_frequency", "frequency"), &LevelDensityGrid::set_path_wobble_frequency);
    ClassDB::bind_method(D_METHOD("get_path_wobble_frequency"), &LevelDensityGrid::get_path_wobble_frequency);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "path_wobble_frequency", PROPERTY_HINT_RANGE, "0.0,1.0,0.01"), "set_path_wobble_frequency", "get_path_wobble_frequency");

    ClassDB::bind_method(D_METHOD("set_connect_from_ground_level", "enabled"), &LevelDensityGrid::set_connect_from_ground_level);
    ClassDB::bind_method(D_METHOD("get_connect_from_ground_level"), &LevelDensityGrid::get_connect_from_ground_level);
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "path_connect_from_ground_level"), "set_connect_from_ground_level", "get_connect_from_ground_level");

    // Smoothing Group
    ADD_GROUP("Smoothing", "smoothing_");
    ClassDB::bind_method(D_METHOD("set_smooth_terrain", "enabled"), &LevelDensityGrid::set_smooth_terrain);
    ClassDB::bind_method(D_METHOD("get_smooth_terrain"), &LevelDensityGrid::get_smooth_terrain);
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "smoothing_enabled"), "set_smooth_terrain", "get_smooth_terrain");
    
    ClassDB::bind_method(D_METHOD("set_smoothing_strength", "strength"), &LevelDensityGrid::set_smoothing_strength);
    ClassDB::bind_method(D_METHOD("get_smoothing_strength"), &LevelDensityGrid::get_smoothing_strength);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "smoothing_strength", PROPERTY_HINT_RANGE, "1,5,1"), "set_smoothing_strength", "get_smoothing_strength");


    ClassDB::bind_method(D_METHOD("_find_ground_position", "start_pos"), &LevelDensityGrid::_find_ground_position);

    ClassDB::bind_method(D_METHOD("generate_liquid_grid"), &LevelDensityGrid::generate_liquid_grid);
    
    ADD_GROUP("Liquid Generation", "liquid_");
    
    // Max Basin Size
    ClassDB::bind_method(D_METHOD("set_max_basin_size", "size"), &LevelDensityGrid::set_max_basin_size);
    ClassDB::bind_method(D_METHOD("get_max_basin_size"), &LevelDensityGrid::get_max_basin_size);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "liquid_max_basin_size", PROPERTY_HINT_RANGE, "1,500,1"), "set_max_basin_size", "get_max_basin_size");

    // Sparsity Cutoff
    ClassDB::bind_method(D_METHOD("set_sparsity_cutoff", "cutoff"), &LevelDensityGrid::set_sparsity_cutoff);
    ClassDB::bind_method(D_METHOD("get_sparsity_cutoff"), &LevelDensityGrid::get_sparsity_cutoff);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "liquid_sparsity_cutoff", PROPERTY_HINT_RANGE, "-1.0,1.0,0.05"), "set_sparsity_cutoff", "get_sparsity_cutoff");

    // Water Height Density
    ClassDB::bind_method(D_METHOD("set_water_height_density", "density"), &LevelDensityGrid::set_water_height_density);
    ClassDB::bind_method(D_METHOD("get_water_height_density"), &LevelDensityGrid::get_water_height_density);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "liquid_water_height_density", PROPERTY_HINT_RANGE, "0.0,1.0,0.01"), "set_water_height_density", "get_water_height_density");

    ClassDB::bind_method(D_METHOD("set_liquid_resolution_multiplier", "multiplier"), &LevelDensityGrid::set_liquid_resolution_multiplier);
    ClassDB::bind_method(D_METHOD("get_liquid_resolution_multiplier"), &LevelDensityGrid::get_liquid_resolution_multiplier);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "liquid_resolution_multiplier", PROPERTY_HINT_RANGE, "1,8,1"), "set_liquid_resolution_multiplier", "get_liquid_resolution_multiplier");

    ClassDB::bind_method(D_METHOD("set_dungeon_path_algorithm", "algo"), &LevelDensityGrid::set_dungeon_path_algorithm);
    ClassDB::bind_method(D_METHOD("get_dungeon_path_algorithm"), &LevelDensityGrid::get_dungeon_path_algorithm);
    
    // Add a dropdown list to the inspector
    ADD_PROPERTY(PropertyInfo(Variant::INT, "dungeon_path_algorithm", PROPERTY_HINT_ENUM, "AStar (Slow/Optimal),Castle (Fast/Recursive)"), "set_dungeon_path_algorithm", "get_dungeon_path_algorithm");
}

// --- Main Generation Function ---

// level_density_grid.cpp

void LevelDensityGrid::generate_from_graph(const TypedArray<Dictionary> &node_data, const TypedArray<Dictionary> &edge_data, float voxel_size, int64_t seed) {
    // 1. Initialize
    initialize_grid(get_grid_size_x(), get_grid_size_y(), get_grid_size_z(), WORLD_SOLID_VALUE);
    resolved_rooms.clear();
    current_carving_zone_id = 0; // Reset active zone
    rng->set_seed(seed);

    // 2. Parse Nodes
    std::vector<ResolvedRoom> processing_rooms;
    std::map<String, int> id_to_index;

    for(int i=0; i < node_data.size(); ++i) {
        Dictionary node = node_data[i];
        ResolvedRoom r;
        r.id = node.get("id", "room");
        r.type = node.get("type", "generic");
        
        Vector3i min_s = node.get("min_size", Vector3i(5,5,5));
        Vector3i max_s = node.get("max_size", Vector3i(10,10,10));
        r.size = Vector3i(
            rng->randi_range(min_s.x, max_s.x),
            rng->randi_range(min_s.y, max_s.y),
            rng->randi_range(min_s.z, max_s.z)
        );

        Vector3i grid_center = get_grid_dimensions() / 2;
        r.position = grid_center + Vector3i(rng->randi_range(-10, 10), 0, rng->randi_range(-10, 10));
        
        processing_rooms.push_back(r);
        id_to_index[r.id] = i;
    }

    // 3. Parse Edges (UPDATED to use ResolvedEdge)
    std::vector<ResolvedEdge> edges;
    for(int i=0; i < edge_data.size(); ++i) {
        Dictionary edge_dict = edge_data[i];
        String from_id = edge_dict.get("from", "");
        String to_id = edge_dict.get("to", "");
        String type = edge_dict.get("type", "corridor"); // Capture the type here!

        if(id_to_index.count(from_id) && id_to_index.count(to_id)) {
            ResolvedEdge e;
            e.from_index = id_to_index[from_id];
            e.to_index = id_to_index[to_id];
            e.type = type;
            edges.push_back(e);
        }
    }

    // 4. Force-Directed Layout
    _resolve_graph_layout(processing_rooms, edges, 100);

    // 5. Carve Rooms
    for(const auto &room : processing_rooms) {
        // Register Zone (NEW)
        int z_id = register_zone_name(room.type);
        current_carving_zone_id = z_id;

        Dictionary room_dict;
        room_dict["start"] = room.position;
        room_dict["end"] = room.position + room.size;
        _create_room(room_dict);
    }

    // 6. Carve Connections
    for(const auto &edge : edges) {
        // Register Zone from the struct (NEW)
        int z_id = register_zone_name(edge.type);
        current_carving_zone_id = z_id;

        ResolvedRoom &rA = processing_rooms[edge.from_index]; // Use .from_index
        ResolvedRoom &rB = processing_rooms[edge.to_index];   // Use .to_index
        
        Vector3 start = rA.center();
        Vector3 end = rB.center();

        if(dungeon_mode) {
             _carve_dungeon_path(start, end);
        } else {
             // Basic bezier placeholder logic
             // In real implementation, pass 'start' and 'end' to your bezier function
             // _carve_bezier_path(start, end); 
             // For now, fallback to simple carve if bezier isn't separate yet:
             _carve_corridor_segment(Vector3i(start), Vector3i(end)); 
        }
    }
    
    // Cleanup
    current_carving_zone_id = 0; 

    // 7. Post-Processing
    if(smooth_terrain) low_pass();
    _calculate_surface_normals();
}

// Updated Layout Solver to handle ResolvedEdge
void LevelDensityGrid::_resolve_graph_layout(std::vector<ResolvedRoom> &rooms, const std::vector<ResolvedEdge> &edges, int iterations) {
    float repulsion = 200.0f;
    float spring_length = 30.0f;
    float spring_strength = 0.5f;
    Vector3 center_target = Vector3(get_grid_size_x()/2, get_grid_size_y()/2, get_grid_size_z()/2);

    for(int k=0; k<iterations; ++k) {
        std::vector<Vector3> forces(rooms.size(), Vector3(0,0,0));

        // Repulsion ... (same as before)
        for(int i=0; i<rooms.size(); ++i) {
            for(int j=i+1; j<rooms.size(); ++j) {
                Vector3 diff = rooms[i].center() - rooms[j].center();
                float dist = diff.length();
                if(dist < 0.1f) dist = 0.1f;
                if (dist < (rooms[i].size.x + rooms[j].size.x)) {
                    Vector3 f = diff.normalized() * (repulsion / dist);
                    forces[i] += f;
                    forces[j] -= f;
                }
            }
            Vector3 to_center = center_target - rooms[i].center();
            forces[i] += to_center * 0.05f; 
        }

        // Springs (Updated accessors)
        for(const auto &edge : edges) {
            int idxA = edge.from_index; // Changed from .first
            int idxB = edge.to_index;   // Changed from .second
            
            Vector3 diff = rooms[idxB].center() - rooms[idxA].center();
            float dist = diff.length();
            float displacement = dist - spring_length;
            Vector3 f = diff.normalized() * (displacement * spring_strength);
            
            forces[idxA] += f;
            forces[idxB] -= f;
        }

        // Apply Forces ... (same as before)
        for(int i=0; i<rooms.size(); ++i) {
            Vector3 new_center = rooms[i].center() + forces[i];
            rooms[i].position = Vector3i(
                new_center.x - rooms[i].size.x / 2,
                Math::clamp((int)new_center.y, 5, get_grid_size_y() - 20),
                new_center.z - rooms[i].size.z / 2
            );
        }
    }
}

void LevelDensityGrid::generate_level_data(const Vector3i &world_grid_dimensions, float voxel_size, int64_t seed) {
    if (world_grid_dimensions.x <= 0 || world_grid_dimensions.y <= 0 || world_grid_dimensions.z <= 0) {
        UtilityFunctions::printerr("LevelDensityGrid.generate_level_data: Invalid World Grid Dimensions ", world_grid_dimensions);
        return;
    }

    UtilityFunctions::print("Generating Level Data for Grid Size: ", world_grid_dimensions, " with seed: ", seed);

    // --- PROFILING START ---
    Time* time = Time::get_singleton();
    uint64_t start_total = time->get_ticks_msec();
    
    uint64_t t_init_grid = 0;
    uint64_t t_rooms_paths = 0;
    uint64_t t_smoothing = 0;
    uint64_t t_normals = 0;

    // 1. Initialize & Seed
    uint64_t t1 = time->get_ticks_usec();
    
    set_noise_seed(seed);
    if (!rng.is_valid()) {
        rng.instantiate();
        UtilityFunctions::printerr("LevelDensityGrid: RNG was invalid, re-initialized.");
    }
    rng->set_seed(seed);

    initialize_grid(world_grid_dimensions.x, world_grid_dimensions.y, world_grid_dimensions.z, WORLD_SOLID_VALUE);
    
    uint64_t t2 = time->get_ticks_usec();
    t_init_grid = (t2 - t1);

    // 2. Generate Rooms and Paths (This includes A* if in dungeon mode)
    _generate_rooms_and_paths(voxel_size);
    
    uint64_t t3 = time->get_ticks_usec();
    t_rooms_paths = (t3 - t2);

    // 3. Optional Smoothing
    if (smooth_terrain && !use_square_brush) {
        UtilityFunctions::print("Smoothing terrain with strength: ", smoothing_strength);
        low_pass();
    }
    
    uint64_t t4 = time->get_ticks_usec();
    t_smoothing = (t4 - t3);

    // 4. Calculate Surface Normals
    _calculate_surface_normals();
    
    uint64_t t5 = time->get_ticks_usec();
    t_normals = (t5 - t4);

    uint64_t end_total = time->get_ticks_msec();

    // --- PRINT REPORT ---
    UtilityFunctions::print("=== Generation Profiling ===");
    UtilityFunctions::print("Total Time:      ", end_total - start_total, " ms");
    UtilityFunctions::print("---------------------------");
    UtilityFunctions::print("Grid Init:       ", t_init_grid / 1000, " ms");
    UtilityFunctions::print("Rooms & Paths:   ", t_rooms_paths / 1000, " ms");
    UtilityFunctions::print("Smoothing:       ", t_smoothing / 1000, " ms");
    UtilityFunctions::print("Surface Normals: ", t_normals / 1000, " ms");
    UtilityFunctions::print("===========================");
}

bool LevelDensityGrid::_is_space_free(const Vector3i& high_res_pos, const PackedFloat32Array& liquid_data, int resolution_mult) {
    // 1. Check Bounds (High Res)
    int hr_gsx = get_grid_size_x() * resolution_mult;
    int hr_gsy = get_grid_size_y() * resolution_mult;
    int hr_gsz = get_grid_size_z() * resolution_mult;

    if (high_res_pos.x < 0 || high_res_pos.x >= hr_gsx ||
        high_res_pos.y < 0 || high_res_pos.y >= hr_gsy ||
        high_res_pos.z < 0 || high_res_pos.z >= hr_gsz) {
        return false;
    }

    // 2. Check Terrain Collision (Map High Res -> Low Res)
    // Integer division handles the mapping (e.g., pos 3 / mult 2 = index 1)
    Vector3i terrain_pos = high_res_pos / resolution_mult;
    
    float terrain_density = get_cell(terrain_pos, 1.0f);
    if (terrain_density > get_surface_threshold()) return false; // Blocked by terrain

    // 3. Check Existing Liquid (High Res)
    // We must manually calculate index because get_index() in this class 
    // is hardcoded for the terrain grid dimensions.
    int64_t idx = (int64_t)high_res_pos.z * hr_gsy * hr_gsx + 
                  (int64_t)high_res_pos.y * hr_gsx + 
                  high_res_pos.x;
                  
    if (idx >= 0 && idx < liquid_data.size()) {
        if (liquid_data[idx] > 0.1f) return false; // Already has water
    }

    return true;
}
// src/level_density_grid.cpp

Ref<DensityGrid> LevelDensityGrid::generate_liquid_grid() {
    Ref<DensityGrid> liquid_grid;
    liquid_grid.instantiate();
    int gsx = get_grid_size_x();
    int gsy = get_grid_size_y();
    int gsz = get_grid_size_z();
    liquid_grid->initialize_grid(gsx, gsy, gsz, 0.0f);

    Ref<FastNoiseLite> distribution_noise;
    distribution_noise.instantiate();
    distribution_noise->set_seed(noise_seed + 999);
    distribution_noise->set_frequency(0.1f); 
    distribution_noise->set_noise_type(FastNoiseLite::TYPE_PERLIN);

    float terrain_thresh = get_surface_threshold();
    
    // USE CLASS MEMBERS HERE INSTEAD OF LOCAL CONSTANTS
    // int max_basin_size = ... (removed)
    // float sparsity_cutoff = ... (removed)

    // Helper: Is this voxel "Solid" enough to hold water?
    auto is_support = [&](Vector3i p) {
        if (p.x < 0 || p.x >= gsx || p.y < 0 || p.y >= gsy || p.z < 0 || p.z >= gsz) return false;
        if (get_cell(p, 1.0f) >= terrain_thresh) return true; 
        if (liquid_grid->get_cell(p, 0.0f) > 0.5f) return true; 
        return false;
    };
    
    auto is_blocked = [&](Vector3i p) {
        if (p.x < 0 || p.x >= gsx || p.y < 0 || p.y >= gsy || p.z < 0 || p.z >= gsz) return true; 
        if (get_cell(p, 1.0f) >= terrain_thresh) return true;
        if (liquid_grid->get_cell(p, 0.0f) > 0.5f) return true;
        return false;
    };

    for (int y = 1; y < gsy - 1; ++y) {
        std::unordered_set<int64_t> layer_visited;
        
        for (int z = 0; z < gsz; ++z) {
            for (int x = 0; x < gsx; ++x) {
                Vector3i pos(x, y, z);
                int64_t idx = liquid_grid->get_index(pos);
                
                if (layer_visited.count(idx)) continue; 
                
                if (is_blocked(pos)) continue; 
                if (!is_support(Vector3i(x, y-1, z))) continue; 
                
                // Use member variable sparsity_cutoff
                if (distribution_noise->get_noise_2d((float)x, (float)z) < sparsity_cutoff) continue;

                std::vector<Vector3i> candidates;
                std::queue<Vector3i> q;
                std::unordered_set<int64_t> current_flood_set; 

                bool leaked = false;
                
                q.push(pos);
                candidates.push_back(pos);
                current_flood_set.insert(idx);
                layer_visited.insert(idx);

                while (!q.empty()) {
                    Vector3i curr = q.front();
                    q.pop();

                    // Use member variable max_basin_size
                    if (candidates.size() > max_basin_size) {
                        leaked = true; 
                        break;
                    }

                    Vector3i nbs[4] = {
                        curr + Vector3i(1,0,0), curr + Vector3i(-1,0,0),
                        curr + Vector3i(0,0,1), curr + Vector3i(0,0,-1)
                    };

                    for (Vector3i nb : nbs) {
                        int64_t nb_idx = liquid_grid->get_index(nb);
                        
                        if (nb.x < 0 || nb.x >= gsx - 1 || nb.z < 0 || nb.z >= gsz - 1 || nb.y < 0 || nb.y >= gsy - 1) {
                            leaked = true;
                            break;
                        }

                        if (is_blocked(nb)) continue;

                        if (!is_support(nb + Vector3i(0, -1, 0))) {
                            leaked = true;
                            break;
                        }

                        if (current_flood_set.find(nb_idx) == current_flood_set.end()) {
                            current_flood_set.insert(nb_idx);
                            layer_visited.insert(nb_idx); 
                            q.push(nb);
                            candidates.push_back(nb);
                        }
                    }
                    if (leaked) break;
                }

                if (!leaked) {
                    for (const Vector3i &p : candidates) {
                        // Use member variable water_height_density
                        liquid_grid->set_cell(p, water_height_density);

                        Vector3i cardinal_dirs[4] = {
                            Vector3i(1, 0, 0), Vector3i(-1, 0, 0),
                            Vector3i(0, 0, 1), Vector3i(0, 0, -1)
                        };

                        for (const Vector3i &dir : cardinal_dirs) {
                            Vector3i nb = p + dir;
                            if (nb.x < 0 || nb.x >= gsx || nb.z < 0 || nb.z >= gsz) continue;

                            if (get_cell(nb, 1.0f) >= terrain_thresh) {
                                // Use member variable water_height_density for bleed
                                liquid_grid->set_cell(nb, water_height_density);
                            }
                        }
                    }
                }
            }
        }
    }

    liquid_grid->set_surface_threshold(0.5f);
    return liquid_grid;
}

// --- Internal Helper Functions ---

void LevelDensityGrid::_calculate_surface_normals() {
    surface_normals.clear();
    int gsx = get_grid_size_x();
    int gsy = get_grid_size_y();
    int gsz = get_grid_size_z();
    float surf_thresh = get_surface_threshold();

    UtilityFunctions::print("Calculating surface normals...");

    for (int x = 0; x < gsx; ++x) {
        for (int y = 0; y < gsy; ++y) {
            for (int z = 0; z < gsz; ++z) {
                Vector3i current_pos(x, y, z);
                if (get_cell(current_pos, WORLD_SOLID_VALUE) > surf_thresh) {
                    bool is_surface = false;
                    Vector3 gradient;

                    // Check X neighbors
                    if (x > 0) {
                        float left = get_cell(Vector3i(x - 1, y, z), WORLD_SOLID_VALUE);
                        if (left < surf_thresh) is_surface = true;
                        gradient.x -= left;
                    }
                    if (x < gsx - 1) {
                        float right = get_cell(Vector3i(x + 1, y, z), WORLD_SOLID_VALUE);
                        if (right < surf_thresh) is_surface = true;
                        gradient.x += right;
                    }
                    // Check Y neighbors
                    if (y > 0) {
                        float down = get_cell(Vector3i(x, y - 1, z), WORLD_SOLID_VALUE);
                        if (down < surf_thresh) is_surface = true;
                        gradient.y -= down;
                    }
                    if (y < gsy - 1) {
                        float up = get_cell(Vector3i(x, y + 1, z), WORLD_SOLID_VALUE);
                        if (up < surf_thresh) is_surface = true;
                        gradient.y += up;
                    }
                    // Check Z neighbors
                    if (z > 0) {
                        float back = get_cell(Vector3i(x, y, z - 1), WORLD_SOLID_VALUE);
                        if (back < surf_thresh) is_surface = true;
                        gradient.z -= back;
                    }
                    if (z < gsz - 1) {
                        float front = get_cell(Vector3i(x, y, z + 1), WORLD_SOLID_VALUE);
                        if (front < surf_thresh) is_surface = true;
                        gradient.z += front;
                    }

                    if (is_surface) {
                        Vector3 normal = gradient.normalized();
                        surface_normals[current_pos] = -normal;
                    }
                }
            }
        }
    }
    UtilityFunctions::print("Finished calculating surface normals. Found ", surface_normals.size(), " surface points.");
}


void LevelDensityGrid::_generate_rooms_and_paths(float voxel_size) {
    if (get_grid_size_x() <= 0) {
        UtilityFunctions::printerr("LevelDensityGrid._generate_rooms_and_paths: Grid not initialized.");
        return;
    }
    if (!rng.is_valid()) {
        UtilityFunctions::printerr("LevelDensityGrid._generate_rooms_and_paths: RNG not valid.");
        return;
    }

    // --- 1. Place Rooms ---
    TypedArray<Dictionary> generated_rooms;
    int placed_room_count = 0;

    UtilityFunctions::print("Attempting to place ", room_count, " rooms.");

    for (int i = 0; i < room_count; ++i) {
        bool room_placed = false;
        for (int try_count = 0; try_count < max_placement_tries; ++try_count) {
            Dictionary new_room = _pick_room();
            if (!new_room.is_empty()) {
                if (!_check_overlap(new_room, generated_rooms)) {
                    _create_room(new_room);
                    generated_rooms.append(new_room);
                    placed_room_count++;
                    room_placed = true;
                    break; // Exit placement tries loop
                }
            }
        }
        if (!room_placed) {
            UtilityFunctions::printerr("Failed to place room ", i + 1, " after ", max_placement_tries, " tries.");
        }
    }

    UtilityFunctions::print("Placed ", generated_rooms.size(), " rooms.");

    if (generated_rooms.size() < 2) {
        UtilityFunctions::printerr("Not enough rooms generated (", generated_rooms.size(), ") to create paths.");
        if (get_grid_size_x() > 0) {
            Vector3i center = get_grid_dimensions() / 2;
            Vector3i spawn_cell = _find_ground_position(center);
            calculated_spawn_position = (Vector3(spawn_cell) + Vector3(0.5f, 0.5f, 0.5f)) * voxel_size;
            calculated_end_position = calculated_spawn_position;
        } else {
            calculated_spawn_position = Vector3(0,0,0);
            calculated_end_position = Vector3(0,0,0);
        }
        return;
    }
    
    // --- 2. Carve Paths Between Rooms ---
    if (dungeon_mode) {
        UtilityFunctions::print("Connecting ", generated_rooms.size(), " rooms with dungeon-style rectangular corridors.");
    } else {
        UtilityFunctions::print("Connecting ", generated_rooms.size(), " rooms with bezier curve paths.");
    }
    
    const float max_x = (float)get_grid_size_x() - 1.0f;
    const float max_y = (float)get_grid_size_y() - 1.0f;
    const float max_z = (float)get_grid_size_z() - 1.0f;

    for (int i = 0; i < generated_rooms.size() - 1; ++i) {
        Dictionary start_room = generated_rooms[i];
        Dictionary end_room = generated_rooms[i + 1];

        Vector3i start_room_start = start_room["start"];
        Vector3i start_room_end = start_room["end"];
        Vector3i end_room_start = end_room["start"];
        Vector3i end_room_end = end_room["end"];

        Vector3 path_start_point;
        Vector3 path_end_point;

        if (connect_from_ground_level) {
            Vector3 start_room_center = Vector3(start_room_start + (start_room_end - start_room_start) / 2);
            Vector3 end_room_center = Vector3(end_room_start + (end_room_end - end_room_start) / 2);
            Vector3 direction_vector = end_room_center - start_room_center;

            int start_y = Math::max(start_room_start.y, 1);
            int start_x, start_z;
            if (abs(direction_vector.x) > abs(direction_vector.z)) {
                start_z = rng->randi_range(start_room_start.z, start_room_end.z - 1);
                start_x = (direction_vector.x > 0) ? start_room_end.x - 1 : start_room_start.x;
            } else {
                start_x = rng->randi_range(start_room_start.x, start_room_end.x - 1);
                start_z = (direction_vector.z > 0) ? start_room_end.z - 1 : start_room_start.z;
            }
            path_start_point = Vector3(start_x, start_y, start_z);

            int end_y = Math::max(end_room_start.y, 1);
            int end_x, end_z;
             if (abs(direction_vector.x) > abs(direction_vector.z)) {
                end_z = rng->randi_range(end_room_start.z, end_room_end.z - 1);
                end_x = (direction_vector.x > 0) ? end_room_start.x : end_room_end.x - 1;
            } else {
                end_x = rng->randi_range(end_room_start.x, end_room_end.x - 1);
                end_z = (direction_vector.z > 0) ? end_room_start.z : end_room_end.z - 1;
            }
            path_end_point = Vector3(end_x, end_y, end_z);

        } else {
            path_start_point = Vector3(start_room_start + (start_room_end - start_room_start) / 2);
            path_end_point = Vector3(end_room_start + (end_room_end - end_room_start) / 2);
        }

        if (dungeon_mode) {
            // [UPDATED] Switch between algorithms
            if (dungeon_path_algorithm == ALGO_CASTLE_RECURSIVE) {
                // New Recursive Logic
                _carve_path_castle(path_start_point, path_end_point);
            } else {
                // Your Existing A* Logic (Legacy)
                _carve_dungeon_path(path_start_point, path_end_point);
            }
        } else {
            std::vector<Vector3> path_nodes;
            path_nodes.push_back(path_start_point);
            int segments = (path_segments > 1) ? path_segments : 1;
            if (segments > 1) {
                Vector3 overall_delta = path_end_point - path_start_point;
                for (int j = 1; j < segments; ++j) {
                    float fraction = static_cast<float>(j) / static_cast<float>(segments);
                    Vector3 base_point = path_start_point + overall_delta * fraction;
                    float offset_mag = (overall_delta.length() / static_cast<float>(segments)) * 0.75f;
                    Vector3 offset(rng->randf_range(-offset_mag, offset_mag), rng->randf_range(-offset_mag, offset_mag), rng->randf_range(-offset_mag, offset_mag));
                    Vector3 intermediate_node = base_point + offset;
                    intermediate_node.x = Math::clamp(intermediate_node.x, 0.0f, max_x);
                    intermediate_node.y = Math::clamp(intermediate_node.y, 0.0f, max_y);
                    intermediate_node.z = Math::clamp(intermediate_node.z, 0.0f, max_z);
                    path_nodes.push_back(intermediate_node);
                }
            }
            path_nodes.push_back(path_end_point);

            for (size_t j = 0; j < path_nodes.size() - 1; ++j) {
                Vector3 p0 = path_nodes[j];
                Vector3 p3 = path_nodes[j + 1];
                Vector3 delta = p3 - p0;
                float offset_magnitude = delta.length() * path_bend_factor;
                Vector3 cp1_offset(rng->randf_range(-offset_magnitude, offset_magnitude), rng->randf_range(-offset_magnitude, offset_magnitude), rng->randf_range(-offset_magnitude, offset_magnitude));
                Vector3 cp2_offset(rng->randf_range(-offset_magnitude, offset_magnitude), rng->randf_range(-offset_magnitude, offset_magnitude), rng->randf_range(-offset_magnitude, offset_magnitude));
                Vector3 p1 = p0 + (delta / 3.0f) + cp1_offset;
                Vector3 p2 = p0 + (delta * 2.0f / 3.0f) + cp2_offset;
                p1.x = Math::clamp(p1.x, 0.0f, max_x); p1.y = Math::clamp(p1.y, 0.0f, max_y); p1.z = Math::clamp(p1.z, 0.0f, max_z);
                p2.x = Math::clamp(p2.x, 0.0f, max_x); p2.y = Math::clamp(p2.y, 0.0f, max_y); p2.z = Math::clamp(p2.z, 0.0f, max_z);
                int num_steps = static_cast<int>(delta.length());
                if (num_steps < 10) num_steps = 10;
                for (int k = 0; k <= num_steps; ++k) {
                    float t = static_cast<float>(k) / static_cast<float>(num_steps);
                    float one_minus_t = 1.0f - t; float t2 = t * t; float t3 = t2 * t; float omt2 = one_minus_t * one_minus_t; float omt3 = omt2 * one_minus_t;
                    Vector3 point_on_curve = (p0 * omt3) + (p1 * 3.0f * omt2 * t) + (p2 * 3.0f * one_minus_t * t2) + (p3 * t3);
                    Vector3 final_point = point_on_curve;
                    if (path_wobble_magnitude > 0.001f) {
                        final_point.x += noise_generator->get_noise_1d(point_on_curve.x * path_wobble_frequency) * path_wobble_magnitude;
                        final_point.y += noise_generator->get_noise_1d(point_on_curve.y * path_wobble_frequency + 1000.0f) * path_wobble_magnitude;
                        final_point.z += noise_generator->get_noise_1d(point_on_curve.z * path_wobble_frequency + 2000.0f) * path_wobble_magnitude;
                    }
                    _mark_brush(Vector3i(final_point), path_brush_min_radius, path_brush_max_radius, WORLD_OPEN_VALUE);
                }
            }
        }
    }

    // --- 3. Calculate Spawn and End Positions ---
    if (generated_rooms.size() > 0) {
        Dictionary first_room = generated_rooms[0];
        Dictionary last_room = generated_rooms[generated_rooms.size() - 1];
        Vector3i start_room_start = first_room["start"];
        Vector3i start_room_end = first_room["end"];
        Vector3i start_spawn_midpoint = start_room_start + (start_room_end - start_room_start) / 2;
        Vector3i spawn_cell = _find_ground_position(start_spawn_midpoint);
        calculated_spawn_position = (Vector3(spawn_cell) + Vector3(0.5f, 0.5f, 0.5f)) * voxel_size;
        Vector3i end_room_start = last_room["start"];
        Vector3i end_room_end = last_room["end"];
        Vector3i end_room_midpoint = end_room_start + (end_room_end - end_room_start) / 2;
        Vector3i end_cell = _find_ground_position(end_room_midpoint);
        calculated_end_position = (Vector3(end_cell) + Vector3(0.5f, 0.5f, 0.5f)) * voxel_size;
    }
    UtilityFunctions::print("Room and Path generation finished.");
}


void LevelDensityGrid::_apply_noise() {
     if (get_grid_size_x() <= 0) {
        UtilityFunctions::printerr("LevelDensityGrid._apply_noise: Grid not initialized.");
        return;
    }
     if (!noise_generator.is_valid()) {
        UtilityFunctions::printerr("LevelDensityGrid._apply_noise: Noise generator is not valid.");
        return;
    }
    UtilityFunctions::print("Applying noise... Scale: ", noise_scale, ", Intensity: ", noise_intensity);
    noise_generator->set_seed(noise_seed);
    noise_generator->set_frequency(noise_frequency);
    float inv_noise_scale = 0.0;
    if (noise_scale > 0.0001) {
        inv_noise_scale = 1.0 / noise_scale;
    } else {
        inv_noise_scale = 1.0;
    }
    int gsx = get_grid_size_x();
    int gsy = get_grid_size_y();
    int gsz = get_grid_size_z();
    float surf_thresh = get_surface_threshold();
    for (int z = 0; z < gsz; ++z) {
        for (int y = 0; y < gsy; ++y) {
            for (int x = 0; x < gsx; ++x) {
                Vector3i pos(x, y, z);
                float noise_val = noise_generator->get_noise_3d(static_cast<double>(x) * inv_noise_scale, static_cast<double>(y) * inv_noise_scale, static_cast<double>(z) * inv_noise_scale);
                float current_cell_value = get_cell(pos, WORLD_SOLID_VALUE);
                float new_density;
                if (current_cell_value < surf_thresh) {
                    float carving_noise = Math::min(noise_val, 0.0f);
                    new_density = current_cell_value + carving_noise * noise_intensity;
                } else {
                    new_density = current_cell_value + noise_val * noise_intensity;
                }
                set_cell(pos, new_density);
            }
        }
    }
    UtilityFunctions::print("Noise application finished.");
}


void LevelDensityGrid::_mark_brush(const Vector3i &center, int radius_lower, int radius_higher, float value) {
    if (!rng.is_valid()) {
        UtilityFunctions::printerr("LevelDensityGrid._mark_brush: RNG not valid.");
        return;
    }
    int r_low = radius_lower;
    int r_high = radius_higher;
    if (r_high < r_low) r_high = r_low;
    if (r_low < 0) r_low = 0;

    int radius = (r_low == r_high) ? r_low : rng->randi_range(r_low, r_high);
    if (radius <= 0) return;

    int radius_sq = radius * radius;
    
    // Grid bounds for clamping
    int gsx = get_grid_size_x();
    int gsy = get_grid_size_y();
    int gsz = get_grid_size_z();

    if(dungeon_mode) {
        for (int i = -radius; i <= radius; ++i) {
            for (int j = 0; j <= radius * 2; ++j) { 
                for (int k = -radius; k <= radius; ++k) {
                    if (use_square_brush || dungeon_mode) {
                        Vector3i current_pos = center + Vector3i(i, j, k);
                        
                        // [FIX] Add explicit boundary checks before setting cell
                        if (current_pos.x > 0 && current_pos.x < gsx - 1 &&
                            current_pos.y > 0 && current_pos.y < gsy - 1 &&
                            current_pos.z > 0 && current_pos.z < gsz - 1) {
                            
                            set_cell(current_pos, value);
                            if (current_carving_zone_id > 0) {
                            // OPTIONAL: Don't overwrite existing zones (keeps Rooms intact if paths cross them)
                                if (get_zone_at(current_pos) == 0) {
                                    set_zone_at(current_pos, current_carving_zone_id);
                                    }
                                }
                            }
                    } else {
                        if (i * i + j * j + k * k < radius_sq) {
                            Vector3i current_pos = center + Vector3i(i, j, k);
                            
                            // [FIX] Add explicit boundary checks
                            if (current_pos.x > 0 && current_pos.x < gsx - 1 &&
                                current_pos.y > 0 && current_pos.y < gsy - 1 &&
                                current_pos.z > 0 && current_pos.z < gsz - 1) {
                                
                                set_cell(current_pos, value);
                                if (current_carving_zone_id > 0) {
                                    // OPTIONAL: Don't overwrite existing zones (keeps Rooms intact if paths cross them)
                                    if (get_zone_at(current_pos) == 0) {
                                        set_zone_at(current_pos, current_carving_zone_id);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else {
        for (int i = -radius; i <= radius; ++i) {
            for (int j = -radius; j <= radius; ++j) {
                for (int k = -radius; k <= radius; ++k) {
                    bool in_shape = false;
                    if (use_square_brush) {
                        in_shape = true; 
                    } else {
                        if (i * i + j * j + k * k < radius_sq) {
                            in_shape = true; 
                        }
                    }

                    if (in_shape) {
                        Vector3i current_pos = center + Vector3i(i, j, k);
                        
                        // [FIX] Replaced simple y>0 check with full boundary check
                        if (current_pos.x > 0 && current_pos.x < gsx - 1 &&
                            current_pos.y > 0 && current_pos.y < gsy - 1 &&
                            current_pos.z > 0 && current_pos.z < gsz - 1) {
                            
                            set_cell(current_pos, value);
                        }
                    }
                }
            }
        }
    }
}


void LevelDensityGrid::_create_room(const Dictionary &room) {
    if (!room.has("start") || !room.has("end")) return;
    
    Vector3i start = room["start"];
    Vector3i end = room["end"];

    // Determine raw bounds
    Vector3i min_corner(Math::min(start.x, end.x), Math::min(start.y, end.y), Math::min(start.z, end.z));
    Vector3i max_corner(Math::max(start.x, end.x), Math::max(start.y, end.y), Math::max(start.z, end.z));

    // Get Grid Dimensions
    int gsx = get_grid_size_x();
    int gsy = get_grid_size_y();
    int gsz = get_grid_size_z();

    // Clamp loops to range [1, size-1] to preserve solid boundary shell
    int start_x = Math::max(min_corner.x, 1);
    int end_x   = Math::min(max_corner.x, gsx - 1);
    
    int start_y = Math::max(min_corner.y, 1); 
    int end_y   = Math::min(max_corner.y, gsy - 1);

    int start_z = Math::max(min_corner.z, 1);
    int end_z   = Math::min(max_corner.z, gsz - 1);

    for (int x = start_x; x < end_x; ++x) {
        for (int y = start_y; y < end_y; ++y) {
            for (int z = start_z; z < end_z; ++z) {
                Vector3i pos(x, y, z);
                
                // 1. Set Density (Carve Air)
                set_cell(pos, WORLD_OPEN_VALUE);
                
                // 2. Set Zone ID (NEW FIX)
                // If we are currently directing a specific room type, stamp it here.
                if (current_carving_zone_id > 0) {
                    set_zone_at(pos, current_carving_zone_id);
                }
            }
        }
    }
}


Dictionary LevelDensityGrid::_pick_room() {
    int gsx = get_grid_size_x();
    int gsy = get_grid_size_y();
    int gsz = get_grid_size_z();
    if (gsx <= 2 || gsy <= 2 || gsz <= 2) return Dictionary();
    if (!rng.is_valid()){
        UtilityFunctions::printerr("LevelDensityGrid._pick_room: RNG not valid.");
        return Dictionary();
    }
    Vector3i valid_min_room_size(Math::clamp(min_room_size.x, 1, gsx - 2), Math::clamp(min_room_size.y, 1, gsy - 2), Math::clamp(min_room_size.z, 1, gsz - 2));
    Vector3i valid_max_room_size(Math::clamp(max_room_size.x, valid_min_room_size.x, gsx - 2), Math::clamp(max_room_size.y, valid_min_room_size.y, gsy - 2), Math::clamp(max_room_size.z, valid_min_room_size.z, gsz - 2));
    if (valid_min_room_size.x > valid_max_room_size.x || valid_min_room_size.y > valid_max_room_size.y || valid_min_room_size.z > valid_max_room_size.z) {
        UtilityFunctions::printerr("LevelDensityGrid._pick_room: MinRoomSize (", valid_min_room_size, ") > MaxRoomSize (", valid_max_room_size, ") after clamping to grid.");
        return Dictionary();
    }
    int size_x = (valid_min_room_size.x == valid_max_room_size.x) ? valid_min_room_size.x : rng->randi_range(valid_min_room_size.x, valid_max_room_size.x);
    int size_y = (valid_min_room_size.y == valid_max_room_size.y) ? valid_min_room_size.y : rng->randi_range(valid_min_room_size.y, valid_max_room_size.y);
    int size_z = (valid_min_room_size.z == valid_max_room_size.z) ? valid_min_room_size.z : rng->randi_range(valid_min_room_size.z, valid_max_room_size.z);
    int max_start_x = gsx - size_x - 1;
    int max_start_y = gsy - size_y - 1;
    int max_start_z = gsz - size_z - 1;
    if (max_start_x < 1 || max_start_y < 1 || max_start_z < 1) {
        return Dictionary();
    }
    Vector3i start(rng->randi_range(1, max_start_x), rng->randi_range(1, max_start_y), rng->randi_range(1, max_start_z));
    Vector3i end = start + Vector3i(size_x, size_y, size_z);
    Dictionary room_data;
    room_data["start"] = start;
    room_data["end"] = end;
    return room_data;
}


bool LevelDensityGrid::_check_overlap(const Dictionary &new_room, const TypedArray<Dictionary> &generated_rooms) {
     if (new_room.is_empty() || !new_room.has("start") || !new_room.has("end")) {
        UtilityFunctions::printerr("_check_overlap: Received invalid new_room.");
        return true;
    }
    Vector3i new_start = new_room["start"];
    Vector3i new_end = new_room["end"];
    for (int i = 0; i < generated_rooms.size(); ++i) {
        Dictionary existing_room = generated_rooms[i];
        if (existing_room.is_empty() || !existing_room.has("start") || !existing_room.has("end")) continue;
        Vector3i existing_start = existing_room["start"];
        Vector3i existing_end = existing_room["end"];
        bool overlap_x = new_start.x < existing_end.x && new_end.x > existing_start.x;
        bool overlap_y = new_start.y < existing_end.y && new_end.y > existing_start.y;
        bool overlap_z = new_start.z < existing_end.z && new_end.z > existing_start.z;
        if (overlap_x && overlap_y && overlap_z) {
            return true;
        }
    }
    return false;
}


Vector3i LevelDensityGrid::_find_ground_position(const Vector3i &start_pos) {
    int gsx = get_grid_size_x();
    int gsy = get_grid_size_y();
    int gsz = get_grid_size_z();
    float surf_thresh = get_surface_threshold();
    Vector3i clamped_start_pos(Math::clamp(start_pos.x, 0, gsx - 1), Math::clamp(start_pos.y, 0, gsy - 1), Math::clamp(start_pos.z, 0, gsz - 1));
    for (int y = clamped_start_pos.y; y >= 0; --y) {
        Vector3i check_pos(clamped_start_pos.x, y, clamped_start_pos.z);
        if (get_cell(check_pos, WORLD_SOLID_VALUE) >= surf_thresh) {
            int ground_surface_y = y + 1;
            ground_surface_y = Math::min(ground_surface_y, gsy - 1);
            return Vector3i(clamped_start_pos.x, ground_surface_y, clamped_start_pos.z);
        }
    }
    return Vector3i(clamped_start_pos.x, 0, clamped_start_pos.z);
}

// --- Helper Struct for A* Pathfinding ---
struct PathNode {
    Vector3i pos;
    Vector3i entered_from; // The direction vector we arrived from
    int g;
    int h;
    bool is_staircase;     // State: Are we currently building a staircase?
    Vector3i parent_pos;   // For path reconstruction

    int f() const { return g + h; }

    // Min-heap priority queue comparator (reversed)
    bool operator>(const PathNode& other) const {
        return f() > other.f();
    }
};

// --- Main Path Carving Function ---
void LevelDensityGrid::_carve_dungeon_path(const Vector3 &start, const Vector3 &end) {
    Time* time = Time::get_singleton();
    uint64_t t_start = time->get_ticks_usec();
    Vector3i start_i = Vector3i(start);
    Vector3i end_i = Vector3i(end);
    
    int gsx = get_grid_size_x();
    int gsy = get_grid_size_y();
    int gsz = get_grid_size_z();
    float surf_thresh = get_surface_threshold();

    // 1. Horizontal Directions (For corridors)
    Vector3i dirs_flat[4] = {
        Vector3i(1, 0, 0), Vector3i(-1, 0, 0),
        Vector3i(0, 0, 1), Vector3i(0, 0, -1)
    };

    // 2. Stair Directions (Combined Horizontal + Vertical)
    // We construct these dynamically during the loop to match the 4 flat directions
    // combined with Up (1) and Down (-1).

    struct DungeonNode {
        Vector3i pos;
        Vector3i entered_from; 
        int g;
        int h; // [NEW] Heuristic cost to target
        
        int straight_dist;      
        int dist_since_stair;   
        int current_stair_len;  
        bool is_in_staircase;   
        
        Vector3i parent_pos;

        // Calculates F-score (Total estimated cost)
        int f() const { return g + h; }

        // Priority Queue puts SMALLEST item at top.
        // So 'greater' operator must return true if THIS > OTHER
        bool operator>(const DungeonNode& other) const { 
            return f() > other.f(); 
        }
    };

    std::priority_queue<DungeonNode, std::vector<DungeonNode>, std::greater<DungeonNode>> open_set;
    
    auto get_idx = [&](const Vector3i& v) -> int64_t { 
        return (int64_t)v.z * gsx * gsy + (int64_t)v.y * gsx + v.x; 
    };
    
    int total_voxels = gsx * gsy * gsz;
    std::vector<int> closed_set_g(total_voxels, 2147483647);

    // --- CONFIGURABLE WEIGHTS ---
    const int COST_MOVE_BASE = 10;
    
    // REDUCED from 500 to 100. Turns are still discouraged, but not "impossible"
    const int COST_TURN_BASE = 100;       
    const int COST_TURN_DECAY = 20;       
    const int COST_TURN_MIN = 10;         

    const int COST_STAIR_START_BASE = 500; // Reduced from 1000
    const int COST_STAIR_START_DECAY = 50;  
    const int COST_STAIR_START_MIN = 50;
    
    const int COST_STAIR_STEP_BASE = 20;    
    const int COST_STAIR_LEN_PENALTY = 10;

    // Reduced from 5000. 
    // If this is too high, it creates "invisible walls" that A* desperately tries to flood-fill around.
    const int COST_ROOM_PENALTY = 500;
    const float HEURISTIC_WEIGHT = 10.0f;

    int start_h = (Math::abs(start_i.x - end_i.x) + 
                   Math::abs(start_i.y - end_i.y) + 
                   Math::abs(start_i.z - end_i.z)) * COST_MOVE_BASE * HEURISTIC_WEIGHT; // <--- Multiply here

    DungeonNode start_node = {
        start_i, Vector3i(0,0,0), 0, start_h, // Pass 'start_h' here
        0, 100, 0, false, start_i
    };
    
    open_set.push(start_node);
    int64_t start_idx = get_idx(start_i);
    closed_set_g[start_idx] = 0;
    std::unordered_map<int64_t, DungeonNode> came_from;

    DungeonNode final_node;
    bool found = false;
    int max_iterations = 2000000; 
    int iter = 0;
    int iterations = 0;

    while (!open_set.empty() && iter < max_iterations) {
        iter++;
        DungeonNode current = open_set.top();
        open_set.pop();

        if (current.pos == end_i) {
            final_node = current;
            found = true;
            break;
        }

        // We check 4 Cardinal directions. 
        // For each cardinal direction, we can go Flat, Up-Stair, or Down-Stair.
        for (const Vector3i& d_flat : dirs_flat) {
            
            // Define the 3 possible moves for this cardinal direction
            struct MoveOption {
                Vector3i dir;
                bool is_stair;
                int y_change;
            };

            MoveOption options[3] = {
                { d_flat, false, 0 },                     // Flat Move
                { d_flat + Vector3i(0, 1, 0), true, 1 },  // Stair Up
                { d_flat + Vector3i(0, -1, 0), true, -1 } // Stair Down
            };

            for (const auto& opt : options) {
                Vector3i next_pos = current.pos + opt.dir;

                // 1. Boundary Check
                if (next_pos.x <= 1 || next_pos.x >= gsx - 2 ||
                    next_pos.y <= 1 || next_pos.y >= gsy - 2 ||
                    next_pos.z <= 1 || next_pos.z >= gsz - 2) {
                    continue;
                }

                // 2. Determine State
                // Note: We compare against d_flat for turns, because even if we go up/down,
                // we are still moving in that cardinal direction.
                bool is_straight = (d_flat == current.entered_from) || (current.entered_from == Vector3i(0,0,0));
                
                // However, if we were in a stair and switch to flat (or vice versa), that's a "state change" cost,
                // effectively acting like a turn or specific penalty.
                
                int move_cost = COST_MOVE_BASE;

                // A. Room Avoidance
                float cell_val = get_cell(next_pos, WORLD_SOLID_VALUE);
                if (cell_val < surf_thresh && next_pos != end_i && next_pos != start_i) {
                    move_cost += COST_ROOM_PENALTY;
                }

                // B. Turn Logic (Horizontal Plane)
                if (!is_straight) {
                    // Turn penalty applies to both flat and stair moves if they change horizontal direction
                    int reduction = current.straight_dist * COST_TURN_DECAY;
                    int turn_cost = Math::max(COST_TURN_MIN, COST_TURN_BASE - reduction);
                    move_cost += turn_cost;
                }

                // C. Stair Logic
                if (opt.is_stair) {
                    if (current.is_in_staircase) {
                        // Continuing a staircase
                        // Check if we are keeping the same vertical direction (Up vs Down)
                        // (Simple check: compare y_change. If we switched Up to Down, huge penalty)
                        Vector3i prev_move = current.pos - current.parent_pos;
                        if (Math::sign(prev_move.y) != opt.y_change) {
                            move_cost += 2000; // Penalty for zigzagging up/down immediately
                        }

                        move_cost += COST_STAIR_STEP_BASE;
                        move_cost += (current.current_stair_len * COST_STAIR_LEN_PENALTY);
                    } else {
                        // STARTING a staircase
                        int reduction = current.dist_since_stair * COST_STAIR_START_DECAY;
                        int start_cost = Math::max(COST_STAIR_START_MIN, COST_STAIR_START_BASE - reduction);
                        move_cost += start_cost;
                    }
                }

                // --- 3. Update Next State ---
                int new_straight_dist = is_straight ? current.straight_dist + 1 : 1;
                int new_dist_since_stair = current.dist_since_stair;
                int new_stair_len = 0;
                
                if (opt.is_stair) {
                    new_stair_len = current.current_stair_len + 1;
                    new_dist_since_stair = 0; // Reset cooldown
                } else {
                    new_stair_len = 0;
                    if (!current.is_in_staircase) {
                         new_dist_since_stair++;
                    }
                }
                
                // Store "entered_from" as the horizontal component (d_flat)
                // so the next iteration knows if we kept going straight horizontally.
                Vector3i new_entered_from = d_flat;

                int new_g = current.g + move_cost;
                int64_t n_idx = get_idx(next_pos);

                if (new_g < closed_set_g[n_idx]) {
                    closed_set_g[n_idx] = new_g;
                    
                    // [NEW] Calculate Heuristic for neighbor
                    int new_h = (Math::abs(next_pos.x - end_i.x) + 
                                 Math::abs(next_pos.y - end_i.y) + 
                                 Math::abs(next_pos.z - end_i.z)) * COST_MOVE_BASE * HEURISTIC_WEIGHT; // <--- Multiply here

                    DungeonNode next_node = {
                        next_pos, new_entered_from, new_g, new_h,
                        new_straight_dist, new_dist_since_stair, new_stair_len, 
                        opt.is_stair, current.pos
                    };
                    came_from[n_idx] = next_node;
                    open_set.push(next_node);
                }
            }
        }
        iterations++;
    }

    uint64_t t_end = time->get_ticks_usec();
    if ((t_end - t_start) > 10000) {
        UtilityFunctions::print("Slow Path Carve: ", (t_end - t_start) / 1000, " ms. Iterations: ", iterations);
    }

    if (found) {
        // Reconstruct
        std::vector<Vector3i> path;
        DungeonNode curr = final_node;
        while (curr.pos != start_i) {
            path.push_back(curr.pos);
            int64_t idx = get_idx(curr.pos);
            if (came_from.find(idx) == came_from.end()) break;
            
            // Parent Lookup
            int64_t parent_idx = get_idx(curr.parent_pos);
            if (came_from.find(parent_idx) != came_from.end()) {
                curr = came_from[parent_idx];
            } else if (curr.parent_pos == start_i) {
                break;
            } else {
                break;
            }
        }
        path.push_back(start_i);

        UtilityFunctions::print("Dungeon Path Generated. Length: ", path.size());

        // Carve
        for (const Vector3i& p : path) {
            // Mark the floor/walking space
            _mark_brush(p, path_brush_min_radius, path_brush_max_radius, WORLD_OPEN_VALUE);
            
            // HEADROOM LOGIC:
            // If this node is part of a staircase (it's not flat relative to neighbors),
            // we must ensure headroom.
            // Since we don't have easy access to neighbors here, we can just aggressively clear upwards
            // for ANY node, or rely on the brush size. 
            // Better: If the path implies verticality, clear 2 blocks up.
            
            // To be safe for navigation, let's clear p + (0,1,0) and p + (0,2,0)
            // explicitly with a smaller brush or exact set_cell to ensure head doesn't clip.
            
            set_cell(p + Vector3i(0, 1, 0), WORLD_OPEN_VALUE);
            set_cell(p + Vector3i(0, 2, 0), WORLD_OPEN_VALUE);
        }
    } else {
        UtilityFunctions::printerr("Dungeon Pathfinding Failed. Falling back.");
        _carve_corridor_segment(start_i, end_i);
    }
}

void LevelDensityGrid::_carve_staircase_v2(const Vector3i &start, const Vector3i &target) {
    Vector3i current = start;
    int y_diff = target.y - start.y;
    int y_dir = Math::sign(y_diff);
    
    // We use a simple 1x1 step logic here
    // For every 1 block we move up/down, we move 1 block forward in a chosen direction
    // to create a 45-degree staircase.
    
    // Choose a direction for the stairs to "unfold" (defaulting to Z if possible)
    Vector3i step_forward(0, 0, 1); 
    
    while (current.y != target.y) {
        // 1. Carve the 'walking' space (the step)
        set_cell(current, WORLD_OPEN_VALUE);
        
        // 2. Clear 2-3 blocks above the step for head clearance
        for (int h = 1; h <= 3; ++h) {
            set_cell(current + Vector3i(0, h, 0), WORLD_OPEN_VALUE);
        }
        
        // 3. Move to the next step (1 forward, 1 up/down)
        current.y += y_dir;
        // current += step_forward; // Optional: include if you want diagonal stairs
    }
    
    // Finish by connecting the end of the stairs to the target room center
    _carve_corridor_segment(current, target);
}

void LevelDensityGrid::_carve_corridor_segment(const Vector3i &from, const Vector3i &to) {
    Vector3i current = from;
    if (current == to) {
        _mark_brush(current, path_brush_min_radius, path_brush_max_radius, WORLD_OPEN_VALUE);
        return;
    }

    // Split diagonal horizontal segments (preserves L-shape logic)
    // Only applies if Y is flat. If Y changes, we rely on the robust loop below.
    if ((from.x != to.x && from.z != to.z) && from.y == to.y) {
        Vector3i mid(to.x, from.y, from.z);
        _carve_corridor_segment(from, mid);
        _carve_corridor_segment(mid, to);
        return;
    }

    _mark_brush(current, path_brush_min_radius, path_brush_max_radius, WORLD_OPEN_VALUE);

    // [FIX] Loop with re-evaluation and safety break
    int sanity_check = 0;
    int max_iter = 10000; // Prevent infinite hangs

    while (current != to && sanity_check < max_iter) {
        
        // Calculate direction for THIS step specifically
        // This handles non-uniform diagonals (e.g. 2 right, 10 up) correctly.
        Vector3i dir(
            (to.x == current.x) ? 0 : static_cast<int>(Math::sign(to.x - current.x)),
            (to.y == current.y) ? 0 : static_cast<int>(Math::sign(to.y - current.y)),
            (to.z == current.z) ? 0 : static_cast<int>(Math::sign(to.z - current.z))
        );

        current += dir;
        _mark_brush(current, path_brush_min_radius, path_brush_max_radius, WORLD_OPEN_VALUE);
        
        sanity_check++;
    }

    if (sanity_check >= max_iter) {
        UtilityFunctions::printerr("LevelDensityGrid: _carve_corridor_segment infinite loop detected from ", from, " to ", to);
    }
}

void LevelDensityGrid::low_pass() {
    int smoothing = smoothing_strength;
    int gsx = get_grid_size_x();
    int gsy = get_grid_size_y();
    int gsz = get_grid_size_z();
    int required_size = smoothing * 2 + 1;
    if (gsx < required_size || gsy < required_size || gsz < required_size) {
        UtilityFunctions::printerr("LevelDensityGrid.low_pass: Grid is too small for the given smoothing strength. Grid: (", get_grid_dimensions(), "), Required: ", required_size);
        return;
    }
    std::vector<float> smoothed_grid;
    size_t grid_total_size = (size_t)gsx * gsy * gsz;
    if (grid_total_size == 0) {
         UtilityFunctions::printerr("LevelDensityGrid.low_pass: Failed to get grid size.");
         return;
    }
    smoothed_grid.resize(grid_total_size);
    for (int z = smoothing; z < gsz - smoothing; ++z) {
        for (int y = smoothing; y < gsy - smoothing; ++y) {
            for (int x = smoothing; x < gsx - smoothing; ++x) {
                Vector3i current_pos(x, y, z);
                float sum = 0.0f;
                int count = 0;
                for (int dz = -smoothing; dz <= smoothing; ++dz) {
                    for (int dy = -smoothing; dy <= smoothing; ++dy) {
                        for (int dx = -smoothing; dx <= smoothing; ++dx) {
                            Vector3i neighbor_pos = current_pos + Vector3i(dx, dy, dz);
                            sum += get_cell(neighbor_pos, WORLD_SOLID_VALUE);
                            count++;
                        }
                    }
                }
                float average = (count > 0) ? sum / count : get_cell(current_pos, WORLD_SOLID_VALUE);
                int index = get_index(current_pos);
                 if (index != -1) {
                    smoothed_grid[index] = average;
                 } else {
                     UtilityFunctions::printerr("LevelDensityGrid.low_pass: Invalid index calculated for ", current_pos);
                 }
            }
        }
    }
    for (int z = smoothing; z < gsz - smoothing; ++z) {
        for (int y = smoothing; y < gsy - smoothing; ++y) {
            for (int x = smoothing; x < gsx - smoothing; ++x) {
                Vector3i pos(x, y, z);
                int index = get_index(pos);
                 if (index != -1) {
                     set_cell(pos, smoothed_grid[index]);
                 }
            }
        }
    }
}

void LevelDensityGrid::_carve_path_castle(const Vector3 &start, const Vector3 &end) {
    // Treat the entire path as one recursive structure. 
    // We pass the full height difference into the recursion.
    _carve_recursive_winding_path(Vector3i(start), Vector3i(end), 4);
}

void LevelDensityGrid::_carve_recursive_winding_path(const Vector3i &start, const Vector3i &end, int depth) {
    Vector3 diff = Vector3(end - start);
    float dist = diff.length();

    // Base Case: Connect directly with a Stepped L-Shape
    if (depth <= 0 || dist < 20.0f) {
        _carve_stepped_L_shape(start, end);
        return;
    }

    // Recursive Step: Midpoint with Height Interpolation
    Vector3i mid = (start + end) / 2;
    
    // [KEY CHANGE] Interpolate height, but snap to integer steps if desired.
    // Simple average distributes height evenly across the path.
    // To make it "occasional", we could add randomness here, but average is safe.
    mid.y = (start.y + end.y) / 2;

    int offset_mag = (int)(dist * 0.3f); 
    if (offset_mag < 2) offset_mag = 2;

    // Displace horizontally
    if (Math::abs(diff.x) > Math::abs(diff.z)) {
        mid.z += rng->randi_range(-offset_mag, offset_mag);
    } else {
        mid.x += rng->randi_range(-offset_mag, offset_mag);
    }
    
    // Clamp X/Z but allow Y to change
    int gsx = get_grid_size_x(); int gsz = get_grid_size_z();
    mid.x = Math::clamp(mid.x, 2, gsx - 3);
    mid.z = Math::clamp(mid.z, 2, gsz - 3);
    
    // Recurse
    _carve_recursive_winding_path(start, mid, depth - 1);
    _carve_recursive_winding_path(mid, end, depth - 1);
}

void LevelDensityGrid::_carve_stepped_L_shape(const Vector3i &start, const Vector3i &end) {
    // We have two legs: Start -> Corner -> End
    // We need to decide which Y level the Corner sits at.
    // It can be at start.y (Leg 1 flat, Leg 2 sloped)
    // OR at end.y (Leg 1 sloped, Leg 2 flat).
    
    bool x_first = rng->randf() > 0.5f;
    bool corner_at_start_height = rng->randf() > 0.5f;

    Vector3i corner_flat; // The X/Z coordinates of the corner
    if (x_first) {
        corner_flat = Vector3i(end.x, 0, start.z);
    } else {
        corner_flat = Vector3i(start.x, 0, end.z);
    }

    Vector3i corner;
    if (corner_at_start_height) {
        corner = Vector3i(corner_flat.x, start.y, corner_flat.z);
        // Leg 1 is Flat
        _carve_variable_height_leg(start, corner);
        // Leg 2 handles the height change
        _carve_variable_height_leg(corner, end);
    } else {
        corner = Vector3i(corner_flat.x, end.y, corner_flat.z);
        // Leg 1 handles the height change
        _carve_variable_height_leg(start, corner);
        // Leg 2 is Flat
        _carve_variable_height_leg(corner, end);
    }
}

// This function carves a straight line that can handle height differences
// by inserting a straight staircase in the middle.
void LevelDensityGrid::_carve_variable_height_leg(const Vector3i &start, const Vector3i &end) {
    int height_diff = end.y - start.y;
    
    // If flat, just use standard carve
    if (height_diff == 0) {
        _carve_corridor_segment(start, end);
        return;
    }
    
    int h_dist = (int)Vector2(start.x - end.x, start.z - end.z).length();
    int stair_len = Math::abs(height_diff); // 1:1 slope
    
    // If the corridor is too short for the stairs, we just carve a direct steep slope
    // (This prevents infinite loops or logic errors, though it might look steep)
    if (stair_len >= h_dist) {
        // Fallback: Just direct carve step-by-step
        _carve_corridor_segment(start, end); 
        return;
    }
    
    // Center the stairs in the corridor
    // Structure: [Flat Part] -> [Stairs] -> [Flat Part]
    int flat_padding = (h_dist - stair_len) / 2;
    
    // Calculate vector direction
    Vector3 dir_vec = Vector3(end - start).normalized();
    // We need precise integer steps for the flat parts
    // It's axis aligned, so dir_vec is clean (e.g. 1,0,0 or 0,0,-1) but Y is tricky.
    // Let's rely on finding intermediate points.
    
    // Find Stair Start point (Flat from Start -> StairStart)
    Vector3i stair_start = start;
    // Advance 'flat_padding' units along X or Z (Y remains start.y)
    if (start.x != end.x) {
        stair_start.x += Math::sign(end.x - start.x) * flat_padding;
    } else {
        stair_start.z += Math::sign(end.z - start.z) * flat_padding;
    }
    
    // Find Stair End point (StairStart -> StairEnd handles height)
    Vector3i stair_end = stair_start;
    if (start.x != end.x) {
        stair_end.x += Math::sign(end.x - start.x) * stair_len;
    } else {
        stair_end.z += Math::sign(end.z - start.z) * stair_len;
    }
    stair_end.y = end.y; // Target height reached here

    // 1. Carve First Flat Section
    _carve_corridor_segment(start, stair_start);
    
    // 2. Carve The Staircase
    // We implement a simple straight stair carver here inline or helper
    // Since it's straight, we iterate stair_len steps.
    Vector3i cursor = stair_start;
    int y_dir = Math::sign(height_diff);
    
    // Determine horizontal step direction
    Vector3i step_dir = Vector3i(0,0,0);
    if (stair_start.x != stair_end.x) step_dir.x = Math::sign(stair_end.x - stair_start.x);
    else step_dir.z = Math::sign(stair_end.z - stair_start.z);

    // Brush settings
    int radius = path_brush_min_radius;
    if (radius < 2) radius = 2;

    for (int i = 0; i < stair_len; ++i) {
        // Carve Step
        _mark_brush(cursor, radius, radius, WORLD_OPEN_VALUE);
        
        // Headroom
        if (radius < 3) _mark_brush(cursor + Vector3i(0, 2, 0), radius, radius, WORLD_OPEN_VALUE);
        
        // Move Up/Down
        cursor.y += y_dir;
        
        // Move Forward
        cursor += step_dir;
    }
    // Ensure the landing is clear
    _mark_brush(cursor, radius, radius, WORLD_OPEN_VALUE);

    // 3. Carve Second Flat Section
    _carve_corridor_segment(stair_end, end);
}

void LevelDensityGrid::set_dungeon_path_algorithm(int p_algo) { dungeon_path_algorithm = p_algo; }
int LevelDensityGrid::get_dungeon_path_algorithm() const { return dungeon_path_algorithm; }


// --- Property Getters/Setters Implementation ---

Vector3 LevelDensityGrid::get_calculated_spawn_position() const { return calculated_spawn_position; }
Vector3 LevelDensityGrid::get_calculated_end_position() const { return calculated_end_position; }
Dictionary LevelDensityGrid::get_surface_normals() const { return surface_normals; }

void LevelDensityGrid::set_noise_scale(float p_scale) { noise_scale = p_scale; }
float LevelDensityGrid::get_noise_scale() const { return noise_scale; }

void LevelDensityGrid::set_noise_intensity(float p_intensity) { noise_intensity = p_intensity; }
float LevelDensityGrid::get_noise_intensity() const { return noise_intensity; }

void LevelDensityGrid::set_noise_seed(int p_seed) {
    noise_seed = p_seed;
    if (noise_generator.is_valid()) {
        noise_generator->set_seed(noise_seed);
    }
}
int LevelDensityGrid::get_noise_seed() const { return noise_seed; }

void LevelDensityGrid::set_noise_frequency(float p_freq) {
    noise_frequency = p_freq;
     if (noise_generator.is_valid()) {
        noise_generator->set_frequency(noise_frequency);
    }
}
float LevelDensityGrid::get_noise_frequency() const { return noise_frequency; }

void LevelDensityGrid::set_room_count(int p_count) { room_count = p_count > 0 ? p_count : 0; }
int LevelDensityGrid::get_room_count() const { return room_count; }

void LevelDensityGrid::set_min_room_size(const Vector3i &p_size) { min_room_size = p_size; }
Vector3i LevelDensityGrid::get_min_room_size() const { return min_room_size; }

void LevelDensityGrid::set_max_room_size(const Vector3i &p_size) { max_room_size = p_size; }
Vector3i LevelDensityGrid::get_max_room_size() const { return max_room_size; }

void LevelDensityGrid::set_max_placement_tries(int p_tries) { max_placement_tries = p_tries > 0 ? p_tries : 1; }
int LevelDensityGrid::get_max_placement_tries() const { return max_placement_tries; }

void LevelDensityGrid::set_path_brush_min_radius(int p_radius) { path_brush_min_radius = p_radius >= 0 ? p_radius : 0; }
int LevelDensityGrid::get_path_brush_min_radius() const { return path_brush_min_radius; }

void LevelDensityGrid::set_path_brush_max_radius(int p_radius) { path_brush_max_radius = p_radius >= 0 ? p_radius : 0; }
int LevelDensityGrid::get_path_brush_max_radius() const { return path_brush_max_radius; }

void LevelDensityGrid::set_use_square_brush(bool p_enabled) { use_square_brush = p_enabled; }
bool LevelDensityGrid::get_use_square_brush() const { return use_square_brush; }

void LevelDensityGrid::set_vertical_movement_cost_multiplier(float p_mult) { vertical_movement_cost_multiplier = p_mult > 0 ? p_mult : 1.0; }
float LevelDensityGrid::get_vertical_movement_cost_multiplier() const { return vertical_movement_cost_multiplier; }

void LevelDensityGrid::set_dungeon_mode(bool p_enabled) { dungeon_mode = p_enabled; }
bool LevelDensityGrid::get_dungeon_mode() const { return dungeon_mode; }

void LevelDensityGrid::set_path_segments(int p_segments) { path_segments = p_segments > 0 ? p_segments : 1; }
int LevelDensityGrid::get_path_segments() const { return path_segments; }

void LevelDensityGrid::set_path_bend_factor(float p_factor) { path_bend_factor = p_factor >= 0.0f ? p_factor : 0.0f; }
float LevelDensityGrid::get_path_bend_factor() const { return path_bend_factor; }

void LevelDensityGrid::set_path_wobble_magnitude(float p_magnitude) { path_wobble_magnitude = p_magnitude >= 0.0f ? p_magnitude : 0.0f; }
float LevelDensityGrid::get_path_wobble_magnitude() const { return path_wobble_magnitude; }

void LevelDensityGrid::set_path_wobble_frequency(float p_frequency) { path_wobble_frequency = p_frequency; }
float LevelDensityGrid::get_path_wobble_frequency() const { return path_wobble_frequency; }

void LevelDensityGrid::set_connect_from_ground_level(bool p_enabled) { connect_from_ground_level = p_enabled; }
bool LevelDensityGrid::get_connect_from_ground_level() const { return connect_from_ground_level; }

void LevelDensityGrid::set_smooth_terrain(bool p_enabled) { smooth_terrain = p_enabled; }
bool LevelDensityGrid::get_smooth_terrain() const { return smooth_terrain; }

void LevelDensityGrid::set_smoothing_strength(int p_strength) { smoothing_strength = p_strength > 0 ? p_strength : 1; }
int LevelDensityGrid::get_smoothing_strength() const { return smoothing_strength; }

void LevelDensityGrid::set_noise_generator(const Ref<FastNoiseLite> &p_noise) {
    noise_generator = p_noise;
    if (noise_generator.is_valid()) {
        noise_generator->set_seed(noise_seed);
        noise_generator->set_frequency(noise_frequency);
    }
}
Ref<FastNoiseLite> LevelDensityGrid::get_noise_generator() const { return noise_generator; }

void LevelDensityGrid::set_max_basin_size(int p_size) { max_basin_size = p_size > 1 ? p_size : 1; }
int LevelDensityGrid::get_max_basin_size() const { return max_basin_size; }

void LevelDensityGrid::set_sparsity_cutoff(float p_cutoff) { sparsity_cutoff = p_cutoff; }
float LevelDensityGrid::get_sparsity_cutoff() const { return sparsity_cutoff; }

void LevelDensityGrid::set_water_height_density(float p_density) { water_height_density = Math::clamp(p_density, 0.0f, 1.0f); }
float LevelDensityGrid::get_water_height_density() const { return water_height_density; }

void LevelDensityGrid::set_liquid_resolution_multiplier(int p_mult) { liquid_resolution_multiplier = (p_mult < 1) ? 1 : p_mult; }
int LevelDensityGrid::get_liquid_resolution_multiplier() const { return liquid_resolution_multiplier; }

String LevelDensityGrid::get_room_type_at(const Vector3 &world_pos, float voxel_size) const {
    // Convert world pos to grid pos
    Vector3i grid_pos = Vector3i(world_pos / voxel_size);
    
    for(const auto &room : resolved_rooms) {
        // Simple AABB check
        if (grid_pos.x >= room.position.x && grid_pos.x < room.position.x + room.size.x &&
            grid_pos.y >= room.position.y && grid_pos.y < room.position.y + room.size.y &&
            grid_pos.z >= room.position.z && grid_pos.z < room.position.z + room.size.z) {
            return room.type;
        }
    }
    return "corridor"; // Default if not in a room
}
} // namespace godot