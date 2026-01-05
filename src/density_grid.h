// src/density_grid.h
#ifndef DENSITY_GRID_H
#define DENSITY_GRID_H

#include <godot_cpp/classes/resource.hpp>
#include <godot_cpp/variant/packed_float32_array.hpp>
#include <godot_cpp/variant/vector3i.hpp>

namespace godot {

class DensityGrid : public Resource {
    GDCLASS(DensityGrid, Resource);

private:
    PackedFloat32Array world_density_grid;
    PackedByteArray world_material_grid;
    PackedInt32Array world_zone_grid;
    std::vector<String> zone_id_map;
    int grid_size_x = 0;
    int grid_size_y = 0;
    int grid_size_z = 0;
    float surface_threshold = 0.0;

protected:
    static void _bind_methods();

public:
    DensityGrid();
    ~DensityGrid();

    // Initialization & Accessors
    void initialize_grid(int size_x, int size_y, int size_z, float default_value = 1.0, int default_material_index = 0);
    bool is_valid_position(const Vector3i &pos) const; // Made const
    int get_index(const Vector3i &pos) const;          // Made const
    bool set_cell(const Vector3i &pos, float value);
    float get_cell(const Vector3i &pos, float default_value = 1.0) const; // Made const
    Vector3i get_grid_dimensions() const; // Made const

    // Properties exposed to Godot
    void set_world_density_grid(const PackedFloat32Array &p_grid);
    PackedFloat32Array get_world_density_grid() const;

    void set_grid_size_x(int p_x);
    int get_grid_size_x() const;

    void set_grid_size_y(int p_y);
    int get_grid_size_y() const;

    void set_grid_size_z(int p_z);
    int get_grid_size_z() const;

    void set_surface_threshold(float p_threshold);
    float get_surface_threshold() const;

    void set_material_id(const Vector3i &pos, int material_index);
    int get_material_id(const Vector3i &pos) const;
    
    // New Property Accessors
    void set_world_material_grid(const PackedByteArray &p_grid);
    PackedByteArray get_world_material_grid() const;

    void set_zone_at(const Vector3i &pos, int zone_id);
    int get_zone_at(const Vector3i &pos) const;
    String get_zone_name_by_id(int zone_id) const;
    int register_zone_name(String name);
};

} // namespace godot

#endif // DENSITY_GRID_H