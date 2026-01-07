// level_density_grid.h
#ifndef LEVEL_DENSITY_GRID_H
#define LEVEL_DENSITY_GRID_H

#include "density_grid.h"
#include "level_gen_data.h"
#include "room_generator.h"
#include "path_carver.h"
#include "liquid_generator.h"

#include <godot_cpp/classes/fast_noise_lite.hpp>
#include <godot_cpp/classes/random_number_generator.hpp>
#include <godot_cpp/variant/vector3.hpp>
#include <godot_cpp/variant/vector3i.hpp>
#include <godot_cpp/variant/typed_array.hpp>
#include <godot_cpp/variant/dictionary.hpp>
#include <godot_cpp/core/object.hpp>

#include "ogt_vox.h"
#include <map>
#include <vector>

namespace godot {

class LevelDensityGrid : public DensityGrid {
    GDCLASS(LevelDensityGrid, DensityGrid);

public:
    float WORLD_OPEN_VALUE = 0.0f;
    float WORLD_SOLID_VALUE = 1.0f;

    // We expose ResolvedRoom/Edge via common header types, 
    // but for GDScript API compatibility we might need to expose them? 
    // Actually structs aren't exposed directly to GDScript usually unless registered.
    // The previous code had them public but they weren't used in API methods (except internal logic).
    
    using ResolvedRoom = godot::ResolvedRoom;
    using ResolvedEdge = godot::ResolvedEdge;

private:
public: // Components Public for easy access internally (or make private and use getters)
    RoomGenerator room_gen;
    PathCarver path_carver;
    LiquidGenerator liquid_gen;

private:
    std::vector<ResolvedRoom> resolved_rooms;
    int current_carving_zone_id = 0;

    // --- Vox Integration ---
    std::map<String, const ogt_vox_scene*> vox_cache;
    
    struct VoxSpawn {
        Vector3 position;
        int palette_index;
    };
    std::vector<VoxSpawn> current_vox_spawns;
    
    // Maps passed from GDScript (cached during generation)
    Dictionary current_vox_spawn_map;
    Dictionary current_vox_material_map;

    const ogt_vox_scene* _load_vox(const String &path);
    void _stamp_vox(const ResolvedRoom &room, const ogt_vox_scene* scene);
    void _get_vox_bounds(const ogt_vox_scene* scene, Vector3i &r_min, Vector3i &r_max);
    void _carve_box(const Vector3i &origin, const Vector3i &size, float value);
    void _clear_vox_cache();


    // --- Configuration Properties (Direct on LevelGrid) ---
    // Noise is distinct from components usually
    float noise_scale = 50.0;
    float noise_intensity = 0.5;
    int noise_seed = 1337;
    float noise_frequency = 0.02;

    bool smooth_terrain = false;
    int smoothing_strength = 1;

    // --- Runtime Data ---
    Vector3 calculated_spawn_position = Vector3(0, 0, 0);
    Vector3 calculated_end_position = Vector3(0, 0, 0);
    Dictionary surface_normals;

    Ref<FastNoiseLite> noise_generator;
    Ref<RandomNumberGenerator> rng;

    // --- Internal Helpers ---
    void _apply_noise(); 
    void _calculate_surface_normals();
    // Helper to bridge the gap or we just expose the component logic via generate_level_data

protected:
    static void _bind_methods();

public:
    LevelDensityGrid();
    ~LevelDensityGrid();

    void generate_level_data(const Vector3i &world_grid_dimensions, float voxel_size, int64_t seed);
    void generate_from_graph(const TypedArray<Dictionary> &node_data, const TypedArray<Dictionary> &edge_data, float voxel_size, int64_t seed, Dictionary settings = Dictionary());
    Vector3 get_calculated_spawn_position() const;
    Vector3 get_calculated_end_position() const;
    Dictionary get_surface_normals() const;

    // Helper for GDScript to mark brush (forwarded to PathCarver logic or reimplemented?)
    // _mark_brush was exposed. It is useful for manual painting. 
    // PathCarver has a private _mark_brush. We can expose a public one here that calls it or duplicates simple logic.
    void _mark_brush(const Vector3i &center, int radius_lower, int radius_higher, float value);
    Vector3i _find_ground_position(const Vector3i &start_pos);

    // --- Getters / Setters (Delegates) ---
    
    // Noise
    void set_noise_scale(float p_scale);
    float get_noise_scale() const;
    void set_noise_intensity(float p_intensity);
    float get_noise_intensity() const;
    void set_noise_seed(int p_seed);
    int get_noise_seed() const;
    void set_noise_frequency(float p_freq);
    float get_noise_frequency() const;
    void set_noise_generator(const Ref<FastNoiseLite> &p_noise);
    Ref<FastNoiseLite> get_noise_generator() const;

    // Room Generator Delegates
    void set_room_count(int p_count);
    int get_room_count() const;
    void set_min_room_size(const Vector3i &p_size);
    Vector3i get_min_room_size() const;
    void set_max_room_size(const Vector3i &p_size);
    Vector3i get_max_room_size() const;
    void set_max_placement_tries(int p_tries);
    int get_max_placement_tries() const;

    // Path Carver Delegates
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
    void set_dungeon_path_algorithm(int p_algo);
    int get_dungeon_path_algorithm() const;
    
    // Bezier
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

    // Smoothing
    void set_smooth_terrain(bool p_enabled);
    bool get_smooth_terrain() const;
    void set_smoothing_strength(int p_strength);
    int get_smoothing_strength() const;
    void low_pass(); // Public API

    // Liquid Generator Delegates
    Ref<DensityGrid> generate_liquid_grid();
    void set_max_basin_size(int p_size);
    int get_max_basin_size() const;
    void set_sparsity_cutoff(float p_cutoff);
    float get_sparsity_cutoff() const;
    void set_water_height_density(float p_density);
    float get_water_height_density() const;
    void set_liquid_resolution_multiplier(int p_mult);
    int get_liquid_resolution_multiplier() const;

    // Vox Integration
    TypedArray<Dictionary> get_vox_spawns() const;
    void set_vox_inverse_density(bool p_enabled);
    bool get_vox_inverse_density() const;
    bool vox_inverse_density = false;

    enum PlacementAlgorithm {
        ALGO_FORCE_DIRECTED = 0,
        ALGO_BSP = 1
    };

    void set_placement_algorithm(PlacementAlgorithm p_algo);
    PlacementAlgorithm get_placement_algorithm() const;

private:
    PlacementAlgorithm placement_algorithm = ALGO_FORCE_DIRECTED;
    void _place_rooms_force_directed(std::vector<ResolvedRoom> &processing_rooms, const std::vector<ResolvedEdge> &edges, const std::vector<int> &node_depths, int max_depth, const TypedArray<Dictionary> &macro_volumes, const PackedVector3Array &guide_points);
    void _place_rooms_bsp(std::vector<ResolvedRoom> &processing_rooms);

    // Level Architecting
    String get_room_type_at(const Vector3 &world_pos, float voxel_size) const;
};

} // namespace godot

VARIANT_ENUM_CAST(godot::LevelDensityGrid::PlacementAlgorithm);

#endif // LEVEL_DENSITY_GRID_H
