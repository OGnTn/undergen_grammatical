#ifndef LEVEL_GEN_DATA_H
#define LEVEL_GEN_DATA_H

#include <godot_cpp/variant/vector3i.hpp>
#include <godot_cpp/variant/vector3.hpp>
#include <godot_cpp/variant/string.hpp>

namespace godot {

struct ResolvedRoom {
    Vector3i position; // Top-left corner in grid
    Vector3i size;
    String type;
    String id;
    
    // Constraints
    bool is_fixed = false;
    Vector3 fixed_position; 
    
    bool has_relative_constraint = false;
    String relative_target_id;
    Vector3 relative_offset;
    
    // Vox Integration
    String vox_path;

    Vector3 center() const { return Vector3(position) + Vector3(size) / 2.0; }
};

struct ResolvedEdge {
    int from_index;
    int to_index;
    String type;
};

} // namespace godot

#endif // LEVEL_GEN_DATA_H
