#ifndef UNDERGEN_LIQUID_FLOOD_NODE_H
#define UNDERGEN_LIQUID_FLOOD_NODE_H

#include "undergen_node.h"
#include "density_grid.h"
#include <godot_cpp/classes/random_number_generator.hpp>
#include <godot_cpp/variant/dictionary.hpp>
#include <godot_cpp/variant/string.hpp>

namespace godot {

class UnderGenLiquidFloodNode : public UnderGenNode {
    GDCLASS(UnderGenLiquidFloodNode, UnderGenNode);

public:
    enum FlowMode {
        FLOW_MODE_IMMEDIATE_POOL = 0,
        FLOW_MODE_MINECRAFT_SIM = 1
    };

private:
    int liquid_material_id = 2;
    Dictionary zone_flood_levels;
    float pool_radius = 5.0f;
    float pool_depth = 3.0f;
    float voxel_size = 1.0f;

    // Minecraft Flow & Basin Carving Properties
    FlowMode flow_mode = FLOW_MODE_IMMEDIATE_POOL;
    int flow_spread_limit = 7;
    float source_spawn_chance = 1.0f;
    int max_sources = 10;
    bool carve_basins = false;
    float basin_radius = 4.0f;
    float basin_depth = 2.0f;
    float basin_carve_value = 0.0f;

    Ref<RandomNumberGenerator> rng;

    // Helper functions for Minecraft Sim & Carving
    void _carve_basin(Ref<DensityGrid> grid, float v_size, const Vector3 &center, float radius, float depth, float surf_thresh, std::vector<Vector3i> &out_basin_cells);
    int _find_distance_to_drop(Ref<DensityGrid> grid, float surf_thresh, const Vector3i &start, int max_dist) const;

protected:
    static void _bind_methods();

public:
    UnderGenLiquidFloodNode();
    virtual ~UnderGenLiquidFloodNode();

    void set_liquid_material_id(int p_id);
    int get_liquid_material_id() const;

    void set_zone_flood_levels(const Dictionary &p_levels);
    Dictionary get_zone_flood_levels() const;

    void set_pool_radius(float p_radius);
    float get_pool_radius() const;

    void set_pool_depth(float p_depth);
    float get_pool_depth() const;

    void set_voxel_size(float p_size);
    float get_voxel_size() const;

    // Flow & Basin Carving Getters/Setters
    void set_flow_mode(FlowMode p_mode);
    FlowMode get_flow_mode() const;

    void set_flow_spread_limit(int p_limit);
    int get_flow_spread_limit() const;

    void set_source_spawn_chance(float p_chance);
    float get_source_spawn_chance() const;

    void set_max_sources(int p_max);
    int get_max_sources() const;

    void set_carve_basins(bool p_carve);
    bool get_carve_basins() const;

    void set_basin_radius(float p_radius);
    float get_basin_radius() const;

    void set_basin_depth(float p_depth);
    float get_basin_depth() const;

    void set_basin_carve_value(float p_value);
    float get_basin_carve_value() const;

    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
};

} // namespace godot

VARIANT_ENUM_CAST(godot::UnderGenLiquidFloodNode::FlowMode);

#endif // UNDERGEN_LIQUID_FLOOD_NODE_H

