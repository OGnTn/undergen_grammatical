#include "undergen_surface_sampler_node.h"
#include "density_grid.h"
#include <godot_cpp/variant/utility_functions.hpp>

namespace godot {

UnderGenSurfaceSamplerNode::UnderGenSurfaceSamplerNode() {}
UnderGenSurfaceSamplerNode::~UnderGenSurfaceSamplerNode() {}

void UnderGenSurfaceSamplerNode::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_surface_type", "type"), &UnderGenSurfaceSamplerNode::set_surface_type);
    ClassDB::bind_method(D_METHOD("get_surface_type"), &UnderGenSurfaceSamplerNode::get_surface_type);
    ClassDB::bind_method(D_METHOD("set_slope_threshold", "threshold"), &UnderGenSurfaceSamplerNode::set_slope_threshold);
    ClassDB::bind_method(D_METHOD("get_slope_threshold"), &UnderGenSurfaceSamplerNode::get_slope_threshold);
    ClassDB::bind_method(D_METHOD("set_voxel_size", "size"), &UnderGenSurfaceSamplerNode::set_voxel_size);
    ClassDB::bind_method(D_METHOD("get_voxel_size"), &UnderGenSurfaceSamplerNode::get_voxel_size);
    ClassDB::bind_method(D_METHOD("set_zone_filter", "filter"), &UnderGenSurfaceSamplerNode::set_zone_filter);
    ClassDB::bind_method(D_METHOD("get_zone_filter"), &UnderGenSurfaceSamplerNode::get_zone_filter);
    ClassDB::bind_method(D_METHOD("set_zone_match_mode", "mode"), &UnderGenSurfaceSamplerNode::set_zone_match_mode);
    ClassDB::bind_method(D_METHOD("get_zone_match_mode"), &UnderGenSurfaceSamplerNode::get_zone_match_mode);

    BIND_ENUM_CONSTANT(FLOOR);
    BIND_ENUM_CONSTANT(CEILING);
    BIND_ENUM_CONSTANT(WALL);
    BIND_ENUM_CONSTANT(ALL);
    BIND_ENUM_CONSTANT(ZONE_MATCH_EXACT);
    BIND_ENUM_CONSTANT(ZONE_MATCH_PREFIX);

    ADD_PROPERTY(PropertyInfo(Variant::INT, "surface_type", PROPERTY_HINT_ENUM, "Floor,Ceiling,Wall,All"),
                 "set_surface_type", "get_surface_type");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "slope_threshold"), "set_slope_threshold", "get_slope_threshold");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "voxel_size"), "set_voxel_size", "get_voxel_size");
    ADD_PROPERTY(PropertyInfo(Variant::STRING, "zone_filter", PROPERTY_HINT_NONE, "Comma-separated zone names. Empty = all."),
                 "set_zone_filter", "get_zone_filter");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "zone_match_mode", PROPERTY_HINT_ENUM, "Exact,Prefix"),
                 "set_zone_match_mode", "get_zone_match_mode");
}

void UnderGenSurfaceSamplerNode::set_surface_type(SurfaceType p_type) { surface_type = p_type; }
UnderGenSurfaceSamplerNode::SurfaceType UnderGenSurfaceSamplerNode::get_surface_type() const { return surface_type; }
void UnderGenSurfaceSamplerNode::set_slope_threshold(float p_threshold) { slope_threshold = p_threshold; }
float UnderGenSurfaceSamplerNode::get_slope_threshold() const { return slope_threshold; }
void UnderGenSurfaceSamplerNode::set_voxel_size(float p_size) { voxel_size = p_size; }
float UnderGenSurfaceSamplerNode::get_voxel_size() const { return voxel_size; }
void UnderGenSurfaceSamplerNode::set_zone_filter(const String &p_filter) { zone_filter = p_filter; }
String UnderGenSurfaceSamplerNode::get_zone_filter() const { return zone_filter; }
void UnderGenSurfaceSamplerNode::set_zone_match_mode(ZoneMatchMode p_mode) { zone_match_mode = p_mode; }
UnderGenSurfaceSamplerNode::ZoneMatchMode UnderGenSurfaceSamplerNode::get_zone_match_mode() const { return zone_match_mode; }

bool UnderGenSurfaceSamplerNode::_zone_matches(const String &point_zone) const {
    if (zone_filter.is_empty()) return true;

    PackedStringArray parts = zone_filter.split(",");
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

void UnderGenSurfaceSamplerNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    Dictionary context = inputs.get(0, Dictionary());
    if (context.is_empty()) return;

    Ref<DensityGrid> grid = context.get("grid", Ref<DensityGrid>());
    if (grid.is_null()) { outputs[0] = context; return; }

    int gsx = grid->get_grid_size_x();
    int gsy = grid->get_grid_size_y();
    int gsz = grid->get_grid_size_z();
    float surf = grid->get_surface_threshold();

    // Build surface normals and generate points
    Ref<UnderGenPointSet> point_set;
    point_set.instantiate();

    static const Vector3i neighbors[6] = {
        Vector3i( 1, 0, 0), Vector3i(-1, 0, 0),
        Vector3i( 0, 1, 0), Vector3i( 0,-1, 0),
        Vector3i( 0, 0, 1), Vector3i( 0, 0,-1)
    };

    for (int z = 0; z < gsz; ++z) {
        for (int y = 0; y < gsy; ++y) {
            for (int x = 0; x < gsx; ++x) {
                Vector3i pos(x, y, z);
                float cell_val = grid->get_cell(pos, 1.0f);
                if (cell_val <= surf) continue; // Must be solid to be a surface

                // Check if any neighbor is air
                bool is_surface = false;
                Vector3 gradient(0, 0, 0);

                for (int n = 0; n < 6; ++n) {
                    Vector3i npos = pos + neighbors[n];
                    float nval = grid->get_cell(npos, 1.0f);
                    if (nval <= surf) is_surface = true;
                    gradient += Vector3(neighbors[n]) * nval;
                }

                if (!is_surface) continue;

                Vector3 normal = -gradient.normalized();

                // Classify surface
                float dot_up = normal.dot(Vector3(0, 1, 0));
                bool is_floor   = dot_up >= slope_threshold;
                bool is_ceiling = dot_up <= -slope_threshold;
                bool is_wall    = !is_floor && !is_ceiling;

                bool passes_filter = false;
                switch (surface_type) {
                    case FLOOR:   passes_filter = is_floor;   break;
                    case CEILING: passes_filter = is_ceiling; break;
                    case WALL:    passes_filter = is_wall;    break;
                    case ALL:     passes_filter = true;       break;
                }

                if (!passes_filter) continue;

                // Zone filter
                if (!_zone_matches(grid->get_zone_name_by_id(grid->get_zone_at(pos)))) continue;

                // Build point attributes
                Dictionary attrs;
                attrs["normal"] = normal;
                attrs["zone"] = grid->get_zone_at(pos);
                attrs["zone_name"] = grid->get_zone_name_by_id(grid->get_zone_at(pos));
                attrs["material"] = grid->get_material_id(pos);
                attrs["slope"] = 1.0f - Math::abs(dot_up); // 0 = flat floor, 1 = vertical wall

                // Offset point slightly above surface in normal direction
                Vector3 world_pos = Vector3(pos) * voxel_size + normal * (voxel_size * 0.5f);
                Transform3D xform;
                xform.origin = world_pos;

                point_set->add_raw_point(xform, 1.0f, attrs);
            }
        }
    }

    outputs[0] = context;           // Port 0: Pass-through context
    outputs[1] = point_set;         // Port 1: Generated PointSet
}

} // namespace godot
