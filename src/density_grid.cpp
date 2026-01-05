// src/density_grid.cpp
#include "density_grid.h"

#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/variant/utility_functions.hpp> // For print/printerr

namespace godot {

DensityGrid::DensityGrid() {
    // Constructor logic if needed
}

DensityGrid::~DensityGrid() {
    // Destructor logic if needed
}

void DensityGrid::_bind_methods() {
    // Methods
    ClassDB::bind_method(D_METHOD("initialize_grid", "size_x", "size_y", "size_z", "default_value", "default_material_index"), &DensityGrid::initialize_grid, DEFVAL(1.0), DEFVAL(0));
    ClassDB::bind_method(D_METHOD("is_valid_position", "pos"), &DensityGrid::is_valid_position);
    ClassDB::bind_method(D_METHOD("get_index", "pos"), &DensityGrid::get_index);
    ClassDB::bind_method(D_METHOD("set_cell", "pos", "value"), &DensityGrid::set_cell);
    ClassDB::bind_method(D_METHOD("get_cell", "pos", "default_value"), &DensityGrid::get_cell, DEFVAL(1.0));
    ClassDB::bind_method(D_METHOD("get_grid_dimensions"), &DensityGrid::get_grid_dimensions);

    // Properties (linked to getters/setters)
    ClassDB::bind_method(D_METHOD("set_world_density_grid", "grid"), &DensityGrid::set_world_density_grid);
    ClassDB::bind_method(D_METHOD("get_world_density_grid"), &DensityGrid::get_world_density_grid);
    ADD_PROPERTY(PropertyInfo(Variant::PACKED_FLOAT32_ARRAY, "world_density_grid"), "set_world_density_grid", "get_world_density_grid");

    ClassDB::bind_method(D_METHOD("set_grid_size_x", "x"), &DensityGrid::set_grid_size_x);
    ClassDB::bind_method(D_METHOD("get_grid_size_x"), &DensityGrid::get_grid_size_x);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "grid_size_x"), "set_grid_size_x", "get_grid_size_x");

    ClassDB::bind_method(D_METHOD("set_grid_size_y", "y"), &DensityGrid::set_grid_size_y);
    ClassDB::bind_method(D_METHOD("get_grid_size_y"), &DensityGrid::get_grid_size_y);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "grid_size_y"), "set_grid_size_y", "get_grid_size_y");

    ClassDB::bind_method(D_METHOD("set_grid_size_z", "z"), &DensityGrid::set_grid_size_z);
    ClassDB::bind_method(D_METHOD("get_grid_size_z"), &DensityGrid::get_grid_size_z);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "grid_size_z"), "set_grid_size_z", "get_grid_size_z");

    ClassDB::bind_method(D_METHOD("set_surface_threshold", "threshold"), &DensityGrid::set_surface_threshold);
    ClassDB::bind_method(D_METHOD("get_surface_threshold"), &DensityGrid::get_surface_threshold);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "surface_threshold"), "set_surface_threshold", "get_surface_threshold");

    ClassDB::bind_method(D_METHOD("set_material_id", "pos", "material_index"), &DensityGrid::set_material_id);
    ClassDB::bind_method(D_METHOD("get_material_id", "pos"), &DensityGrid::get_material_id);

    // 3. Bind new property
    ClassDB::bind_method(D_METHOD("set_world_material_grid", "grid"), &DensityGrid::set_world_material_grid);
    ClassDB::bind_method(D_METHOD("get_world_material_grid"), &DensityGrid::get_world_material_grid);
    ADD_PROPERTY(PropertyInfo(Variant::PACKED_BYTE_ARRAY, "world_material_grid"), "set_world_material_grid", "get_world_material_grid");

    ClassDB::bind_method(D_METHOD("get_zone_name_by_id", "id"), &DensityGrid::get_zone_name_by_id);
    ClassDB::bind_method(D_METHOD("get_zone_at", "pos"), &DensityGrid::get_zone_at);
}


void DensityGrid::initialize_grid(int size_x, int size_y, int size_z, float default_value, int default_material_index) {
    if (size_x <= 0 || size_y <= 0 || size_z <= 0) {
        UtilityFunctions::printerr("DensityGrid.initialize_grid: Invalid dimensions provided (", size_x, ", ", size_y, ", ", size_z, ").");
        grid_size_x = 0;
        grid_size_y = 0;
        grid_size_z = 0;
        world_density_grid.clear(); // Use clear instead of resize(0)
        world_material_grid.clear();
        return;
    }

    grid_size_x = size_x;
    grid_size_y = size_y;
    grid_size_z = size_z;

    int64_t total_grid_points = (int64_t)grid_size_x * grid_size_y * grid_size_z; // Use int64_t for potentially large sizes
    world_density_grid.resize(total_grid_points);
    world_density_grid.fill(default_value); // fill exists for PackedFloat32Array

    world_material_grid.resize(total_grid_points);
    world_material_grid.fill(default_material_index);

    world_zone_grid.resize(total_grid_points);
    world_zone_grid.fill(0); // 0 = Unclaimed/Solid
    zone_id_map.clear();
    zone_id_map.push_back("undefined");

    UtilityFunctions::print("Initialized density grid of size ", grid_size_x, " x ", grid_size_y, " x ", grid_size_z, " = ", total_grid_points, " points with value ", default_value);
}

bool DensityGrid::is_valid_position(const Vector3i &pos) const {
    return pos.x >= 0 && pos.x < grid_size_x &&
           pos.y >= 0 && pos.y < grid_size_y &&
           pos.z >= 0 && pos.z < grid_size_z;
}

int DensityGrid::get_index(const Vector3i &pos) const {
    if (!is_valid_position(pos)) {
        return -1; // Indicate invalid index
    }
    // Calculate 1D index - Use int64_t for intermediate calculation to prevent overflow
    int64_t index = (int64_t)pos.x + (int64_t)grid_size_x * (pos.y + (int64_t)grid_size_y * pos.z);
    // Check if the calculated index fits within a 32-bit int range if needed, though PackedArray size implies it should
    if (index < 0 || index >= world_density_grid.size()) {
         // This case should technically not happen if is_valid_position is correct and grid is initialized
        // UtilityFunctions::printerr("DensityGrid.get_index: Calculated index out of bounds (", index, ") for pos ", pos);
        return -1;
    }
    return (int)index;
}


bool DensityGrid::set_cell(const Vector3i &pos, float value) {
    int index = get_index(pos);
    if (index != -1) {
        // Ensure index is within the valid range (double check)
        if (index >= 0 && index < world_density_grid.size()) {
             world_density_grid[index] = value;
             return true;
        } else {
            // UtilityFunctions::printerr("DensityGrid.set_cell: Index ", index, " out of bounds [0, ", world_density_grid.size(), ") despite valid position ", pos);
             return false; // Should not happen if get_index is correct
        }
    }
    // Optionally log warning for OOB write attempts
    // UtilityFunctions::printerr("DensityGrid.set_cell: Attempted write out of bounds at ", pos, " or grid not initialized.");
    return false;
}


float DensityGrid::get_cell(const Vector3i &pos, float default_value) const {
    int index = get_index(pos);
    if (index != -1) {
         // Ensure index is within the valid range (double check)
        if (index >= 0 && index < world_density_grid.size()) {
            return world_density_grid[index];
        } else {
           // UtilityFunctions::printerr("DensityGrid.get_cell: Index ", index, " out of bounds [0, ", world_density_grid.size(), ") despite valid position ", pos);
            return default_value; // Should not happen if get_index is correct
        }
    }
    // Return default value if out of bounds or grid not initialized
    return default_value;
}

Vector3i DensityGrid::get_grid_dimensions() const {
    return Vector3i(grid_size_x, grid_size_y, grid_size_z);
}

// Property Getters/Setters
void DensityGrid::set_world_density_grid(const PackedFloat32Array &p_grid) {
    world_density_grid = p_grid;
    // Note: Grid dimensions are NOT automatically updated here.
    // The user should call initialize_grid or set dimensions manually
    // if they set the density array directly.
}

PackedFloat32Array DensityGrid::get_world_density_grid() const {
    return world_density_grid;
}

void DensityGrid::set_material_id(const Vector3i &pos, int material_index) {
    int index = get_index(pos);
    if (index != -1 && index < world_material_grid.size()) {
        // Clamp to Byte range (0-255)
        if (material_index < 0) material_index = 0;
        if (material_index > 255) material_index = 255;
        world_material_grid[index] = (uint8_t)material_index;
    }
}

int DensityGrid::get_material_id(const Vector3i &pos) const {
    int index = get_index(pos);
    if (index != -1 && index < world_material_grid.size()) {
        return (int)world_material_grid[index];
    }
    return 0; // Default fallback
}

void DensityGrid::set_world_material_grid(const PackedByteArray &p_grid) {
    world_material_grid = p_grid;
}

PackedByteArray DensityGrid::get_world_material_grid() const {
    return world_material_grid;
}

void DensityGrid::set_grid_size_x(int p_x) { grid_size_x = p_x > 0 ? p_x : 0; }
int DensityGrid::get_grid_size_x() const { return grid_size_x; }

void DensityGrid::set_grid_size_y(int p_y) { grid_size_y = p_y > 0 ? p_y : 0; }
int DensityGrid::get_grid_size_y() const { return grid_size_y; }

void DensityGrid::set_grid_size_z(int p_z) { grid_size_z = p_z > 0 ? p_z : 0; }
int DensityGrid::get_grid_size_z() const { return grid_size_z; }

void DensityGrid::set_surface_threshold(float p_threshold) { surface_threshold = p_threshold; }
float DensityGrid::get_surface_threshold() const { return surface_threshold; }

void DensityGrid::set_zone_at(const Vector3i &pos, int zone_id) {
    int index = get_index(pos);
    if (index != -1 && index < world_zone_grid.size()) {
        world_zone_grid[index] = zone_id;
    }
}

int DensityGrid::get_zone_at(const Vector3i &pos) const {
    int index = get_index(pos);
    if (index != -1 && index < world_zone_grid.size()) {
        return world_zone_grid[index];
    }
    return 0;
}

int DensityGrid::register_zone_name(String name) {
    zone_id_map.push_back(name);
    return zone_id_map.size() - 1;
}

String DensityGrid::get_zone_name_by_id(int zone_id) const {
    if (zone_id > 0 && zone_id < zone_id_map.size()) {
        return zone_id_map[zone_id];
    }
    return "";
}


} // namespace godot