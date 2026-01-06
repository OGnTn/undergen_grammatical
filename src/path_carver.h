#ifndef PATH_CARVER_H
#define PATH_CARVER_H

#include "level_gen_data.h"
#include "density_grid.h"
#include <godot_cpp/classes/fast_noise_lite.hpp>
#include <godot_cpp/classes/random_number_generator.hpp>
#include <godot_cpp/variant/dictionary.hpp>
#include <vector>
#include <queue>
#include <map>

namespace godot {

class PathCarver {
public:
    enum DungeonPathAlgorithm {
        ALGO_ASTAR = 0,
        ALGO_CASTLE_RECURSIVE = 1
    };

    // Configuration
    int path_brush_min_radius = 2;
    int path_brush_max_radius = 4;
    bool use_square_brush = false;
    float vertical_movement_cost_multiplier = 2.0;
    bool dungeon_mode = false; 
    
    // Bezier specific
    int path_segments = 1;
    float path_bend_factor = 0.4f;
    float path_wobble_magnitude = 0.0f;
    float path_wobble_frequency = 0.2f;
    bool connect_from_ground_level = false;

    // AStar specific
    int dungeon_path_algorithm = ALGO_ASTAR;

    // State
    int current_carving_zone_id = 0;

    void connect_rooms(DensityGrid* grid, RandomNumberGenerator* rng, Ref<FastNoiseLite> wobble_noise, const TypedArray<Dictionary>& generated_rooms, float voxel_size, Vector3& out_spawn, Vector3& out_end);
    
    void create_paths_from_edges(DensityGrid* grid, RandomNumberGenerator* rng, Ref<FastNoiseLite> wobble_noise, const std::vector<ResolvedRoom>& rooms, const std::vector<ResolvedEdge>& edges);

private:
    float WORLD_OPEN_VALUE = 0.0f;
    float WORLD_SOLID_VALUE = 1.0f;

    struct DungeonNode {
        Vector3i pos;
        Vector3i arrived_from_dir; // e.g. (1,0,0) if we walked +X to get here
        int g;
        int h; // Heuristic
        int straight_dist;     // Distance moved in same direction
        int dist_since_stair;  // Distance moved since last stair step
        int current_stair_len; // Length of current stair sequence
        bool is_in_staircase;  // Are we currently on a stair (slope)?
        Vector3i parent_pos;
        
        bool operator>(const DungeonNode& other) const {
            return (g + h) > (other.g + other.h); // Min-heap priority
        }
    };

    void _mark_brush(DensityGrid* grid, RandomNumberGenerator* rng, const Vector3i &center, int radius_lower, int radius_higher, float value);
    Vector3i _find_ground_position(DensityGrid* grid, const Vector3i &start_pos);
    
    // Dungeon Carving
    void _carve_dungeon_path(DensityGrid* grid, const Vector3 &start, const Vector3 &end);
    void _carve_bezier_path(DensityGrid* grid, RandomNumberGenerator* rng, Ref<FastNoiseLite> wobble_noise, const Vector3 &start, const Vector3 &end);
    void _carve_corridor_segment(DensityGrid* grid, RandomNumberGenerator* rng, const Vector3i &from, const Vector3i &to);
    
    // Recursive / Castle Algo
    void _carve_path_castle(DensityGrid* grid, RandomNumberGenerator* rng, const Vector3 &start, const Vector3 &end);
    void _carve_recursive_winding_path(DensityGrid* grid, RandomNumberGenerator* rng, const Vector3i &start, const Vector3i &end, int depth);
    void _carve_stepped_L_shape(DensityGrid* grid, RandomNumberGenerator* rng, const Vector3i &start, const Vector3i &end);
    void _carve_variable_height_leg(DensityGrid* grid, RandomNumberGenerator* rng, const Vector3i &start, const Vector3i &end);
    
    // Bezier
    // (Logic integrated into connect_rooms directly in original code)
};

} // namespace godot

#endif // PATH_CARVER_H
