#ifndef UNDERGEN_SURFACE_SAMPLER_NODE_H
#define UNDERGEN_SURFACE_SAMPLER_NODE_H

#include "undergen_node.h"
#include "undergen_point_set.h"
#include <godot_cpp/variant/vector3.hpp>
#include <godot_cpp/variant/vector3i.hpp>
#include <godot_cpp/variant/string.hpp>
#include <godot_cpp/variant/packed_string_array.hpp>

namespace godot {

class UnderGenSurfaceSamplerNode : public UnderGenNode {
    GDCLASS(UnderGenSurfaceSamplerNode, UnderGenNode);

public:
    enum SurfaceType {
        FLOOR = 0,  // Normal pointing up (Y+)
        CEILING = 1, // Normal pointing down (Y-)
        WALL = 2,   // Normal mostly horizontal
        ALL = 3
    };

    enum ZoneMatchMode {
        ZONE_MATCH_EXACT = 0,
        ZONE_MATCH_PREFIX = 1
    };

private:
    SurfaceType surface_type = FLOOR;
    float slope_threshold = 0.6f; // Dot product threshold to classify floor vs wall
    float voxel_size = 1.0f;

    // Zone filter — comma-separated. Empty = sample all zones.
    String zone_filter = "";
    ZoneMatchMode zone_match_mode = ZONE_MATCH_EXACT;

    bool _zone_matches(const String &point_zone) const;

protected:
    static void _bind_methods();

public:
    UnderGenSurfaceSamplerNode();
    virtual ~UnderGenSurfaceSamplerNode();

    void set_surface_type(SurfaceType p_type);
    SurfaceType get_surface_type() const;
    void set_slope_threshold(float p_threshold);
    float get_slope_threshold() const;
    void set_voxel_size(float p_size);
    float get_voxel_size() const;

    void set_zone_filter(const String &p_filter);
    String get_zone_filter() const;
    void set_zone_match_mode(ZoneMatchMode p_mode);
    ZoneMatchMode get_zone_match_mode() const;

    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
};

} // namespace godot

VARIANT_ENUM_CAST(godot::UnderGenSurfaceSamplerNode::SurfaceType);
VARIANT_ENUM_CAST(godot::UnderGenSurfaceSamplerNode::ZoneMatchMode);

#endif // UNDERGEN_SURFACE_SAMPLER_NODE_H
