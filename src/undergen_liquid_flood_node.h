#ifndef UNDERGEN_LIQUID_FLOOD_NODE_H
#define UNDERGEN_LIQUID_FLOOD_NODE_H

#include "undergen_node.h"
#include <godot_cpp/variant/dictionary.hpp>
#include <godot_cpp/variant/string.hpp>

namespace godot {

class UnderGenLiquidFloodNode : public UnderGenNode {
    GDCLASS(UnderGenLiquidFloodNode, UnderGenNode);

private:
    int liquid_material_id = 2;
    Dictionary zone_flood_levels;
    float pool_radius = 5.0f;
    float pool_depth = 3.0f;
    float voxel_size = 1.0f;

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

    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
};

} // namespace godot

#endif // UNDERGEN_LIQUID_FLOOD_NODE_H
