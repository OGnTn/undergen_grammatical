#ifndef ROOM_GENERATOR_H
#define ROOM_GENERATOR_H

#include "level_gen_data.h"
#include "density_grid.h"
#include <godot_cpp/classes/random_number_generator.hpp>
#include <godot_cpp/variant/dictionary.hpp>

namespace godot {

class RoomGenerator {
public:
    // Configuration
    int room_count = 10;
    Vector3i min_room_size = Vector3i(5, 5, 3);
    Vector3i max_room_size = Vector3i(15, 15, 8);
    int max_placement_tries = 5000;
    float voxel_size = 1.0f; // For converting grid pos to world pos

    // Dependencies (Passed during generation usually, but we can store them if needed or pass to method)
    // For this design, we will pass them to the generate method.

    // Logic
    TypedArray<Dictionary> generate_rooms(DensityGrid* grid, RandomNumberGenerator* rng, int& out_placed_count);
    
    // Helper to architect rooms from graph data
    void create_rooms_from_data(DensityGrid* grid, const std::vector<ResolvedRoom>& rooms, int& current_zone_id);

private:
    float WORLD_OPEN_VALUE = 0.0f;

    Dictionary _pick_room(const Vector3i& grid_size, RandomNumberGenerator* rng);
    bool _check_overlap(const Dictionary &new_room, const TypedArray<Dictionary> &generated_rooms);
    void _create_room(DensityGrid* grid, const Dictionary &room, int zone_id);
};

} // namespace godot

#endif // ROOM_GENERATOR_H
