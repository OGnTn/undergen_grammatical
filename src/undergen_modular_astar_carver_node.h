#ifndef UNDERGEN_MODULAR_ASTAR_CARVER_NODE_H
#define UNDERGEN_MODULAR_ASTAR_CARVER_NODE_H

#include "undergen_node.h"
#include "path_carver.h"
#include "vox_spawn_entry.h"
#include "vox_material_entry.h"
#include <godot_cpp/classes/fast_noise_lite.hpp>
#include <godot_cpp/classes/random_number_generator.hpp>
#include <godot_cpp/variant/dictionary.hpp>
#include <godot_cpp/variant/typed_array.hpp>
#include <map>

// Forward declare ogt_vox structures to keep header clean
struct ogt_vox_scene;

namespace godot {

class UnderGenModularAStarCarverNode : public UnderGenNode {
    GDCLASS(UnderGenModularAStarCarverNode, UnderGenNode);

private:
    // Pathfinding Configuration
    int path_brush_min_radius = 2;
    int path_brush_max_radius = 4;
    bool use_square_brush = false;
    float vertical_movement_cost_multiplier = 2.0f;
    bool connect_from_ground_level = false;

    // VOX Tile Paths
    String doorway_vox_path = "";
    String straight_vox_path = "";
    String turn_vox_path = "";
    String stairs_vox_path = "";
    String t_junction_vox_path = "";

    // VOX Stamp Customization
    int connection_palette_index = 0;
    bool vox_inverse_density = false;
    bool exclude_from_smoothing = false;
    bool exclude_from_warping = false;

    TypedArray<VoxSpawnEntry> vox_spawn_entries;
    TypedArray<VoxMaterialEntry> vox_material_entries;
    Dictionary vox_spawn_map;
    Dictionary vox_material_map;

    Ref<RandomNumberGenerator> rng;
    Ref<FastNoiseLite> wobble_noise;

    // VOX Cache for performance
    std::map<String, const ogt_vox_scene*> vox_cache;

    const ogt_vox_scene* _load_vox(const String &path);
    void _clear_vox_cache();
    void _rebuild_maps_from_entries();
    void _on_entries_changed();

    void _stamp_vox_tile(DensityGrid* grid, const Vector3i &grid_pos, const ogt_vox_scene* scene, int zone_id, int rotation_idx, std::vector<Dictionary> &out_spawns, Array &out_rooms);
    void _stamp_vox_tile_fallback(DensityGrid* grid, const Vector3i& center, int zone_id);

protected:
    static void _bind_methods();

public:
    UnderGenModularAStarCarverNode();
    virtual ~UnderGenModularAStarCarverNode();

    // Setters/Getters for Pathfinding
    void set_path_brush_min_radius(int p_radius);
    int get_path_brush_min_radius() const;
    void set_path_brush_max_radius(int p_radius);
    int get_path_brush_max_radius() const;
    void set_use_square_brush(bool p_enabled);
    bool get_use_square_brush() const;
    void set_vertical_movement_cost_multiplier(float p_mult);
    float get_vertical_movement_cost_multiplier() const;
    void set_connect_from_ground_level(bool p_enabled);
    bool get_connect_from_ground_level() const;

    // Setters/Getters for VOX Tile Paths
    void set_doorway_vox_path(const String &p_path);
    String get_doorway_vox_path() const;
    void set_straight_vox_path(const String &p_path);
    String get_straight_vox_path() const;
    void set_turn_vox_path(const String &p_path);
    String get_turn_vox_path() const;
    void set_stairs_vox_path(const String &p_path);
    String get_stairs_vox_path() const;
    void set_t_junction_vox_path(const String &p_path);
    String get_t_junction_vox_path() const;

    // Setters/Getters for VOX Stamp configurations
    void set_connection_palette_index(int p_index);
    int get_connection_palette_index() const;
    void set_vox_inverse_density(bool p_enabled);
    bool get_vox_inverse_density() const;
    void set_exclude_from_smoothing(bool p_exclude);
    bool get_exclude_from_smoothing() const;
    void set_exclude_from_warping(bool p_exclude);
    bool get_exclude_from_warping() const;

    void set_vox_spawn_map(const Dictionary &p_map);
    Dictionary get_vox_spawn_map() const;
    void set_vox_material_map(const Dictionary &p_map);
    Dictionary get_vox_material_map() const;

    void set_vox_spawn_entries(const TypedArray<VoxSpawnEntry> &p_entries);
    TypedArray<VoxSpawnEntry> get_vox_spawn_entries() const;
    void set_vox_material_entries(const TypedArray<VoxMaterialEntry> &p_entries);
    TypedArray<VoxMaterialEntry> get_vox_material_entries() const;

    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
};

} // namespace godot

#endif // UNDERGEN_MODULAR_ASTAR_CARVER_NODE_H
