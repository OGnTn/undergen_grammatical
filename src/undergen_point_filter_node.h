#ifndef UNDERGEN_POINT_FILTER_NODE_H
#define UNDERGEN_POINT_FILTER_NODE_H

#include "undergen_node.h"
#include "undergen_point_set.h"
#include <godot_cpp/variant/string.hpp>

namespace godot {

// Filters a PointSet by attribute conditions, density, zone, slope, and minimum spacing.
class UnderGenPointFilterNode : public UnderGenNode {
    GDCLASS(UnderGenPointFilterNode, UnderGenNode);

public:
    enum ZoneMatchMode {
        ZONE_MATCH_EXACT = 0,         // Point passes if zone matches any name in list
        ZONE_MATCH_PREFIX = 1,        // Point passes if zone starts with any prefix in list
        ZONE_MATCH_EXCLUDE = 2        // Point passes if zone matches NONE of the names in list
    };

private:
    // Zone filter — comma-separated list of zone names. Empty = accept all.
    String required_zone_name = "";
    ZoneMatchMode zone_match_mode = ZONE_MATCH_EXACT;
    // Slope filter (0=flat, 1=vertical)
    float min_slope = 0.0f;
    float max_slope = 1.0f;
    // Density threshold - points below this density value are discarded
    float min_density = 0.0f;
    // Minimum spacing between kept points (in world units). 0 = no spacing check
    float min_spacing = 0.0f;

    bool _zone_matches(const String &point_zone) const;

protected:
    static void _bind_methods();

public:
    UnderGenPointFilterNode();
    virtual ~UnderGenPointFilterNode();

    void set_required_zone_name(const String &p_zone);
    String get_required_zone_name() const;
    void set_zone_match_mode(ZoneMatchMode p_mode);
    ZoneMatchMode get_zone_match_mode() const;
    void set_min_slope(float p_slope);
    float get_min_slope() const;
    void set_max_slope(float p_slope);
    float get_max_slope() const;
    void set_min_density(float p_density);
    float get_min_density() const;
    void set_min_spacing(float p_spacing);
    float get_min_spacing() const;

    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
};

} // namespace godot

VARIANT_ENUM_CAST(godot::UnderGenPointFilterNode::ZoneMatchMode);

#endif // UNDERGEN_POINT_FILTER_NODE_H
