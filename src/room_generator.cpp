#include "room_generator.h"
#include <godot_cpp/variant/utility_functions.hpp>
#include <godot_cpp/core/math.hpp>

namespace godot {

TypedArray<Dictionary> RoomGenerator::generate_rooms(DensityGrid* grid, RandomNumberGenerator* rng, int& out_placed_count) {
    TypedArray<Dictionary> generated_rooms;
    out_placed_count = 0;

    if (!grid) return generated_rooms;
    
    Vector3i grid_dim = grid->get_grid_dimensions();

    UtilityFunctions::print("RoomGenerator: Attempting to place ", room_count, " rooms.");

    for (int i = 0; i < room_count; ++i) {
        bool room_placed = false;
        for (int try_count = 0; try_count < max_placement_tries; ++try_count) {
            Dictionary new_room = _pick_room(grid_dim, rng);
            if (!new_room.is_empty()) {
                if (!_check_overlap(new_room, generated_rooms)) {
                    _create_room(grid, new_room, 0); // Default zone 0? Logic suggests standard rooms don't need zone IDs potentially, or we need to pass it?
                    // Original code: _create_room(new_room) used member current_carving_zone_id implicitly (likely 0) or explicit from graph.
                    // For random generation, usually 0 unless we want to tag them. Let's keep 0 for now.
                    
                    generated_rooms.append(new_room);
                    out_placed_count++;
                    room_placed = true;
                    break; 
                }
            }
        }
        if (!room_placed) {
            UtilityFunctions::printerr("Failed to place room ", i + 1, " after ", max_placement_tries, " tries.");
        }
    }

    UtilityFunctions::print("Placed ", generated_rooms.size(), " rooms.");
    return generated_rooms;
}

void RoomGenerator::create_rooms_from_data(DensityGrid* grid, const std::vector<ResolvedRoom>& rooms, int& current_zone_id) {
    for(const auto &room : rooms) {
        // Here we might want to register zones via the grid (if it handled that) or just pass the ID.
        // The original code called grid->register_zone_name(room.type).
        // Since RoomGenerator doesn't own zone registration, we assume the caller handles zone name reg
        // and sets current_zone_id. 
        // WAIT: The loop was:
        // int z_id = register_zone_name(room.type);
        // current_carving_zone_id = z_id;
        // ...
        
        // This means we need the Grid's help to get the ID.
        int z_id = grid->register_zone_name(room.type);
        
        Dictionary room_dict;
        room_dict["start"] = room.position;
        room_dict["end"] = room.position + room.size;
        
        _create_room(grid, room_dict, z_id);
    }
}

Dictionary RoomGenerator::_pick_room(const Vector3i& grid_size, RandomNumberGenerator* rng) {
    int gsx = grid_size.x;
    int gsy = grid_size.y;
    int gsz = grid_size.z;
    if (gsx <= 2 || gsy <= 2 || gsz <= 2) return Dictionary();
    
    // Clamp room sizes to grid
    Vector3i valid_min_room_size(
        Math::clamp(min_room_size.x, 1, gsx - 2), 
        Math::clamp(min_room_size.y, 1, gsy - 2), 
        Math::clamp(min_room_size.z, 1, gsz - 2)
    );
    Vector3i valid_max_room_size(
        Math::clamp(max_room_size.x, valid_min_room_size.x, gsx - 2), 
        Math::clamp(max_room_size.y, valid_min_room_size.y, gsy - 2), 
        Math::clamp(max_room_size.z, valid_min_room_size.z, gsz - 2)
    );

    if (valid_min_room_size.x > valid_max_room_size.x || 
        valid_min_room_size.y > valid_max_room_size.y || 
        valid_min_room_size.z > valid_max_room_size.z) {
        UtilityFunctions::printerr("RoomGenerator._pick_room: MinRoomSize > MaxRoomSize after clamping.");
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

    Vector3i start(
        rng->randi_range(1, max_start_x), 
        rng->randi_range(1, max_start_y), 
        rng->randi_range(1, max_start_z)
    );
    Vector3i end = start + Vector3i(size_x, size_y, size_z);
    
    Dictionary room_data;
    room_data["start"] = start;
    room_data["end"] = end;
    return room_data;
}

bool RoomGenerator::_check_overlap(const Dictionary &new_room, const TypedArray<Dictionary> &generated_rooms) {
     if (new_room.is_empty() || !new_room.has("start") || !new_room.has("end")) {
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

void RoomGenerator::_create_room(DensityGrid* grid, const Dictionary &room, int zone_id) {
    if (!room.has("start") || !room.has("end")) return;
    
    Vector3i start = room["start"];
    Vector3i end = room["end"];

    Vector3i min_corner(Math::min(start.x, end.x), Math::min(start.y, end.y), Math::min(start.z, end.z));
    Vector3i max_corner(Math::max(start.x, end.x), Math::max(start.y, end.y), Math::max(start.z, end.z));

    int gsx = grid->get_grid_size_x();
    int gsy = grid->get_grid_size_y();
    int gsz = grid->get_grid_size_z();

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
                
                grid->set_cell(pos, WORLD_OPEN_VALUE);
                
                if (zone_id > 0) {
                    grid->set_zone_at(pos, zone_id);
                }
            }
        }
    }
}

} // namespace godot
