// level_density_grid.h
#ifndef LEVEL_DENSITY_GRID_H
#define LEVEL_DENSITY_GRID_H

#include "density_grid.h"
#include <godot_cpp/classes/fast_noise_lite.hpp>
#include <godot_cpp/classes/random_number_generator.hpp>
#include <godot_cpp/variant/vector3.hpp>
#include <godot_cpp/variant/vector3i.hpp>
#include <godot_cpp/variant/typed_array.hpp>
#include <godot_cpp/variant/dictionary.hpp>
#include <godot_cpp/core/object.hpp>

namespace godot {

class LevelDensityGrid : public DensityGrid {
    GDCLASS(LevelDensityGrid, DensityGrid);

public:
    float WORLD_OPEN_VALUE = 0.0f;
    float WORLD_SOLID_VALUE = 1.0f;
    enum DungeonPathAlgorithm {
        ALGO_ASTAR = 0,
        ALGO_CASTLE_RECURSIVE = 1
    };

private:
    // --- Configuration Properties ---
    float noise_scale = 50.0;
    float noise_intensity = 0.5;
    int noise_seed = 1337;
    float noise_frequency = 0.02;

    int room_count = 10;
    Vector3i min_room_size = Vector3i(5, 5, 3);
    Vector3i max_room_size = Vector3i(15, 15, 8);
    int max_placement_tries = 5000;

    int path_brush_min_radius = 2;
    int path_brush_max_radius = 4;
    bool use_square_brush = false;
    float vertical_movement_cost_multiplier = 2.0;
    bool dungeon_mode = false; 
    
    int path_segments = 1;
    float path_bend_factor = 0.4f;
    float path_wobble_magnitude = 0.0f;
    float path_wobble_frequency = 0.2f;
    bool connect_from_ground_level = false;

    bool smooth_terrain = false;
    int smoothing_strength = 1;

    int max_basin_size = 40;       // Default from our testing
    float sparsity_cutoff = 0.2f;  // Default noise threshold
    float water_height_density = 0.5f; // Default density for 80% height

    int liquid_resolution_multiplier = 2; // Default to 2x resolution
    int dungeon_path_algorithm = ALGO_ASTAR;

    // Update helper signature to accept the multiplier for terrain lookups
    bool _is_space_free(const Vector3i& high_res_pos, const PackedFloat32Array& liquid_data, int resolution_mult);

    // --- Runtime Data ---
    Vector3 calculated_spawn_position = Vector3(0, 0, 0);
    Vector3 calculated_end_position = Vector3(0, 0, 0);
    Dictionary surface_normals;

    Ref<FastNoiseLite> noise_generator;
    Ref<RandomNumberGenerator> rng;

    // --- Internal Helper Functions ---
    void _generate_rooms_and_paths(float voxel_size);
    void _apply_noise();

    void _create_room(const Dictionary &room);
    Dictionary _pick_room();
    bool _check_overlap(const Dictionary &new_room, const TypedArray<Dictionary> &generated_rooms);
    void _calculate_surface_normals();

    // --- Refined Dungeon Path Helpers ---
    void _carve_dungeon_path(const Vector3 &start, const Vector3 &end);
    void _carve_corridor_segment(const Vector3i &from, const Vector3i &to);
    // V2 implementation for specific stepped stair generation
    void _carve_staircase_v2(const Vector3i &start, const Vector3i &target);
    void _carve_path_castle(const Vector3 &start, const Vector3 &end);
    void _carve_recursive_winding_path(const Vector3i &start, const Vector3i &end, int depth);
    void _carve_stepped_L_shape(const Vector3i &start, const Vector3i &end); // Replaces _carve_L_shape
    void _carve_variable_height_leg(const Vector3i &start, const Vector3i &end); // Replaces vertical_staircase
    
protected:
    static void _bind_methods();

public:
    LevelDensityGrid();
    ~LevelDensityGrid();

    void generate_level_data(const Vector3i &world_grid_dimensions, float voxel_size, int64_t seed);

    Vector3 get_calculated_spawn_position() const;
    Vector3 get_calculated_end_position() const;
    Dictionary get_surface_normals() const;

    void _mark_brush(const Vector3i &center, int radius_lower, int radius_higher, float value);
    Vector3i _find_ground_position(const Vector3i &start_pos);
    void set_noise_scale(float p_scale);
    float get_noise_scale() const;
    void set_noise_intensity(float p_intensity);
    float get_noise_intensity() const;
    void set_noise_seed(int p_seed);
    int get_noise_seed() const;
    void set_noise_frequency(float p_freq);
    float get_noise_frequency() const;

    void set_room_count(int p_count);
    int get_room_count() const;
    void set_min_room_size(const Vector3i &p_size);
    Vector3i get_min_room_size() const;
    void set_max_room_size(const Vector3i &p_size);
    Vector3i get_max_room_size() const;
    void set_max_placement_tries(int p_tries);
    int get_max_placement_tries() const;
    void low_pass();

    void set_path_brush_min_radius(int p_radius);
    int get_path_brush_min_radius() const;
    void set_path_brush_max_radius(int p_radius);
    int get_path_brush_max_radius() const;
    void set_use_square_brush(bool p_enabled);
    bool get_use_square_brush() const;
    void set_vertical_movement_cost_multiplier(float p_mult);
    float get_vertical_movement_cost_multiplier() const;
    void set_dungeon_mode(bool p_enabled);
    bool get_dungeon_mode() const;
    
    void set_path_segments(int p_segments);
    int get_path_segments() const;
    void set_path_bend_factor(float p_factor);
    float get_path_bend_factor() const;
    void set_path_wobble_magnitude(float p_magnitude);
    float get_path_wobble_magnitude() const;
    void set_path_wobble_frequency(float p_frequency);
    float get_path_wobble_frequency() const;

    void set_connect_from_ground_level(bool p_enabled);
    bool get_connect_from_ground_level() const;

    void set_smooth_terrain(bool p_enabled);
    bool get_smooth_terrain() const;
    void set_smoothing_strength(int p_strength);
    int get_smoothing_strength() const;

    void set_noise_generator(const Ref<FastNoiseLite> &p_noise);
    Ref<FastNoiseLite> get_noise_generator() const;
    Ref<DensityGrid> generate_liquid_grid();

    void set_max_basin_size(int p_size);
    int get_max_basin_size() const;

    void set_sparsity_cutoff(float p_cutoff);
    float get_sparsity_cutoff() const;

    void set_water_height_density(float p_density);
    float get_water_height_density() const;

    void set_liquid_resolution_multiplier(int p_mult);
    int get_liquid_resolution_multiplier() const;

    void set_dungeon_path_algorithm(int p_algo);
    int get_dungeon_path_algorithm() const;
};

} // namespace godot

#endif // LEVEL_DENSITY_GRID_H