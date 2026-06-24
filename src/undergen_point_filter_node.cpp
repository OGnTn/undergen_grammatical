#include "undergen_point_filter_node.h"
#include <godot_cpp/core/math.hpp>
#include <godot_cpp/variant/packed_string_array.hpp>
#include <vector>

namespace godot {

UnderGenPointFilterNode::UnderGenPointFilterNode() {}
UnderGenPointFilterNode::~UnderGenPointFilterNode() {}

void UnderGenPointFilterNode::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_required_zone_name", "zone"), &UnderGenPointFilterNode::set_required_zone_name);
    ClassDB::bind_method(D_METHOD("get_required_zone_name"), &UnderGenPointFilterNode::get_required_zone_name);
    ClassDB::bind_method(D_METHOD("set_zone_match_mode", "mode"), &UnderGenPointFilterNode::set_zone_match_mode);
    ClassDB::bind_method(D_METHOD("get_zone_match_mode"), &UnderGenPointFilterNode::get_zone_match_mode);
    ClassDB::bind_method(D_METHOD("set_min_slope", "slope"), &UnderGenPointFilterNode::set_min_slope);
    ClassDB::bind_method(D_METHOD("get_min_slope"), &UnderGenPointFilterNode::get_min_slope);
    ClassDB::bind_method(D_METHOD("set_max_slope", "slope"), &UnderGenPointFilterNode::set_max_slope);
    ClassDB::bind_method(D_METHOD("get_max_slope"), &UnderGenPointFilterNode::get_max_slope);
    ClassDB::bind_method(D_METHOD("set_min_density", "density"), &UnderGenPointFilterNode::set_min_density);
    ClassDB::bind_method(D_METHOD("get_min_density"), &UnderGenPointFilterNode::get_min_density);
    ClassDB::bind_method(D_METHOD("set_min_spacing", "spacing"), &UnderGenPointFilterNode::set_min_spacing);
    ClassDB::bind_method(D_METHOD("get_min_spacing"), &UnderGenPointFilterNode::get_min_spacing);

    BIND_ENUM_CONSTANT(ZONE_MATCH_EXACT);
    BIND_ENUM_CONSTANT(ZONE_MATCH_PREFIX);

    ADD_PROPERTY(PropertyInfo(Variant::STRING, "required_zone_name", PROPERTY_HINT_NONE, "Comma-separated zone names"), "set_required_zone_name", "get_required_zone_name");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "zone_match_mode", PROPERTY_HINT_ENUM, "Exact,Prefix"), "set_zone_match_mode", "get_zone_match_mode");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "min_slope"), "set_min_slope", "get_min_slope");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "max_slope"), "set_max_slope", "get_max_slope");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "min_density"), "set_min_density", "get_min_density");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "min_spacing"), "set_min_spacing", "get_min_spacing");
}

void UnderGenPointFilterNode::set_required_zone_name(const String &p_zone) { required_zone_name = p_zone; }
String UnderGenPointFilterNode::get_required_zone_name() const { return required_zone_name; }
void UnderGenPointFilterNode::set_zone_match_mode(ZoneMatchMode p_mode) { zone_match_mode = p_mode; }
UnderGenPointFilterNode::ZoneMatchMode UnderGenPointFilterNode::get_zone_match_mode() const { return zone_match_mode; }
void UnderGenPointFilterNode::set_min_slope(float p_slope) { min_slope = p_slope; }
float UnderGenPointFilterNode::get_min_slope() const { return min_slope; }
void UnderGenPointFilterNode::set_max_slope(float p_slope) { max_slope = p_slope; }
float UnderGenPointFilterNode::get_max_slope() const { return max_slope; }
void UnderGenPointFilterNode::set_min_density(float p_density) { min_density = p_density; }
float UnderGenPointFilterNode::get_min_density() const { return min_density; }
void UnderGenPointFilterNode::set_min_spacing(float p_spacing) { min_spacing = p_spacing; }
float UnderGenPointFilterNode::get_min_spacing() const { return min_spacing; }

bool UnderGenPointFilterNode::_zone_matches(const String &point_zone) const {
    if (required_zone_name.is_empty()) return true;

    // Split comma-separated filter list
    PackedStringArray parts = required_zone_name.split(",");
    for (int i = 0; i < parts.size(); ++i) {
        String filter = parts[i].strip_edges();
        if (filter.is_empty()) continue;

        switch (zone_match_mode) {
            case ZONE_MATCH_EXACT:
                if (point_zone == filter) return true;
                break;
            case ZONE_MATCH_PREFIX:
                if (point_zone.begins_with(filter)) return true;
                break;
        }
    }
    return false;
}

void UnderGenPointFilterNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    // Port 0: Incoming PointSet
    Ref<UnderGenPointSet> in_points = inputs.get(0, Ref<UnderGenPointSet>());
    if (in_points.is_null()) return;

    Ref<UnderGenPointSet> out_points;
    out_points.instantiate();

    float spacing_sq = min_spacing * min_spacing;
    std::vector<Vector3> kept_positions; // For spacing checks

    const auto& raw = in_points->get_points_raw();
    for (const auto& p : raw) {
        // Density filter
        if (p.density < min_density) continue;

        const Dictionary& attrs = p.attributes;

        // Zone filter (supports comma-separated list + prefix matching)
        if (!_zone_matches(attrs.get("zone_name", ""))) continue;

        // Slope filter
        float slope = attrs.get("slope", 0.0f);
        if (slope < min_slope || slope > max_slope) continue;

        // Spacing filter (Poisson-disk-like rejection)
        if (min_spacing > 0.0f) {
            Vector3 pos = p.transform.origin;
            bool too_close = false;
            for (const Vector3& kp : kept_positions) {
                if (pos.distance_squared_to(kp) < spacing_sq) {
                    too_close = true;
                    break;
                }
            }
            if (too_close) continue;
            kept_positions.push_back(pos);
        }

        out_points->add_point(p.transform, p.density, p.bounds, p.color, p.attributes);
    }

    outputs[0] = out_points; // Port 0: Filtered PointSet
}

} // namespace godot
