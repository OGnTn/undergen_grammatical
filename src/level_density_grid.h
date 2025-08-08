// level_density_grid.h
#ifndef LEVEL_DENSITY_GRID_H
#define LEVEL_DENSITY_GRID_H

#include "density_grid.h" // Include the base class
#include <godot_cpp/classes/fast_noise_lite.hpp>
#include <godot_cpp/classes/a_star3d.hpp>
#include <godot_cpp/classes/random_number_generator.hpp> // For RNG
#include <godot_cpp/variant/vector3.hpp>
#include <godot_cpp/variant/vector3i.hpp>
#include <godot_cpp/variant/typed_array.hpp> // For Array[Dictionary]
#include <godot_cpp/variant/dictionary.hpp>  // For Dictionary
#include <godot_cpp/core/object.hpp> // For Ref<>

namespace godot {

class LevelDensityGrid : public DensityGrid {
    GDCLASS(LevelDensityGrid, DensityGrid);

public: // Constants accessible from outside if needed
    float WORLD_OPEN_VALUE = 0.0f;
    float WORLD_SOLID_VALUE = 1.0f;



private:
    // --- Configuration Properties ---
    // Noise
    float noise_scale = 50.0;
    float noise_intensity = 0.5;
    int noise_seed = 1337;
    float noise_frequency = 0.02;

    // Rooms
    int room_count = 10;
    Vector3i min_room_size = Vector3i(5, 5, 3);
    Vector3i max_room_size = Vector3i(15, 15, 8);
    int max_placement_tries = 5000;

    // Paths
    int path_brush_min_radius = 2;
    int path_brush_max_radius = 4;
    bool use_square_brush = false;              // NEW: Toggle for square brush shape.
    float vertical_movement_cost_multiplier = 2.0; // Note: AStar3D doesn't directly use this per-edge
    bool dungeon_mode = false;                  // Toggle for dungeon-style corridors.
    
    // Parameters for Bezier curve generation
    int path_segments = 1;                  // Number of bezier segments per connection. >1 creates "windy" paths.
    float path_bend_factor = 0.4f;          // How far control points can stray. Controls "bendyness".
    float path_wobble_magnitude = 0.0f;     // Magnitude of high-frequency noise on path. Controls "wobble".
    float path_wobble_frequency = 0.2f;     // Frequency of high-frequency noise on path.
    bool connect_from_ground_level = false; // Toggle to connect rooms from walls at ground level.

    // Smoothing
    bool smooth_terrain = false;
    int smoothing_strength = 1;


    // --- Runtime Data ---
    Vector3 calculated_spawn_position = Vector3(0, 0, 0);
    Vector3 calculated_end_position = Vector3(0, 0, 0);
    Dictionary surface_normals;

    // --- Noise Generator ---
    Ref<FastNoiseLite> noise_generator;

    // --- RNG ---
    Ref<RandomNumberGenerator> rng;


    // --- Internal Helper Functions ---
    void _generate_rooms_and_paths(float voxel_size);
    void _apply_noise();
    void _mark_brush(const Vector3i &center, int radius_lower, int radius_higher, float value);
    void _create_room(const Dictionary &room);
    Dictionary _pick_room(); // Returns Dictionary (empty if failed)
    bool _check_overlap(const Dictionary &new_room, const TypedArray<Dictionary> &generated_rooms);
    void _calculate_surface_normals();
    // Helpers for Dungeon Mode
    void _carve_dungeon_path(const Vector3 &start, const Vector3 &end);
    void _carve_line(const Vector3i &from, const Vector3i &to);
    void _carve_staircase(const Vector3i &bottom, const Vector3i &top, const Vector3i &prev_corner);
    

protected:
    static void _bind_methods();

public:
    LevelDensityGrid();
    ~LevelDensityGrid();

    // --- Main Generation Function ---
    void generate_level_data(const Vector3i &world_grid_dimensions, float voxel_size, int64_t seed);

    // --- Getters for Runtime Data ---
    Vector3 get_calculated_spawn_position() const;
	Vector3 get_calculated_end_position() const;
    Dictionary get_surface_normals() const;

    // --- Property Getters/Setters ---
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
    
    // Bezier Path Getters/Setters
    void set_path_segments(int p_segments);
    int get_path_segments() const;
    void set_path_bend_factor(float p_factor);
    float get_path_bend_factor() const;
    void set_path_wobble_magnitude(float p_magnitude);
    float get_path_wobble_magnitude() const;
    void set_path_wobble_frequency(float p_frequency);
    float get_path_wobble_frequency() const;

    // Ground Connection Getter/Setter
    void set_connect_from_ground_level(bool p_enabled);
    bool get_connect_from_ground_level() const;

    // Smoothing Getters/Setters
    void set_smooth_terrain(bool p_enabled);
    bool get_smooth_terrain() const;
    void set_smoothing_strength(int p_strength);
    int get_smoothing_strength() const;

    void set_noise_generator(const Ref<FastNoiseLite> &p_noise);
    Ref<FastNoiseLite> get_noise_generator() const;
};

} // namespace godot

#endif // LEVEL_DENSITY_GRID_H