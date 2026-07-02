#include "undergen_surface_sampler_node.h"
#include "density_grid.h"
#include "grid_parallel.h"
#include <godot_cpp/variant/utility_functions.hpp>
#include <vector>

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
    int64_t total_size = grid->get_total_cell_count();
    if (total_size <= 0 ||
        grid->get_density_data().size() < total_size ||
        grid->get_zone_data().size() < total_size ||
        grid->get_material_data().size() < total_size) {
        outputs[0] = context;
        return;
    }

    float surf = grid->get_surface_threshold();
    const PackedFloat32Array &density_array = grid->get_density_data();
    const PackedInt32Array &zone_array = grid->get_zone_data();
    const PackedByteArray &material_array = grid->get_material_data();
    const float *density_data = density_array.ptr();
    const int32_t *zone_data = zone_array.ptr();
    const uint8_t *material_data = material_array.ptr();

    // Build surface normals and generate points
    Ref<UnderGenPointSet> point_set;
    point_set.instantiate();

    struct Neighbor {
        int dx;
        int dy;
        int dz;
    };
    static const Neighbor neighbors[6] = {
        { 1, 0, 0}, {-1, 0, 0},
        { 0, 1, 0}, { 0,-1, 0},
        { 0, 0, 1}, { 0, 0,-1}
    };

    std::vector<String> zone_names(grid->get_zone_count());
    for (int i = 0; i < (int)zone_names.size(); ++i) {
        zone_names[i] = grid->get_zone_name_by_id(i);
    }

    std::vector<String> filters;
    if (!zone_filter.is_empty()) {
        PackedStringArray parts = zone_filter.split(",");
        for (int i = 0; i < parts.size(); ++i) {
            String filter = parts[i].strip_edges();
            if (!filter.is_empty()) {
                filters.push_back(filter);
            }
        }
    }

    auto zone_matches_fast = [&](const String &point_zone) -> bool {
        if (filters.empty()) return true;

        for (const String &filter : filters) {
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
    };

    struct SampledPoint {
        Vector3 world_pos;
        Vector3 normal;
        int zone = 0;
        String zone_name;
        int material = 0;
        float slope = 0.0f;
    };

    int worker_count = grid_parallel_worker_count(gsz, total_size);
    std::vector<std::vector<SampledPoint>> worker_points(worker_count);

    parallel_for_z(gsz, total_size, [&](int worker_index, int z_begin, int z_end) {
    std::vector<SampledPoint> &local_points = worker_points[worker_index];
    for (int z = z_begin; z < z_end; ++z) {
        int slice_offset = z * gsy * gsx;
        for (int y = 0; y < gsy; ++y) {
            int row_offset = slice_offset + y * gsx;
            for (int x = 0; x < gsx; ++x) {
                int idx = row_offset + x;
                float cell_val = density_data[idx];
                if (cell_val <= surf) continue; // Must be solid to be a surface

                // Check if any neighbor is air
                bool is_surface = false;
                Vector3 gradient(0, 0, 0);

                for (int n = 0; n < 6; ++n) {
                    int nx = x + neighbors[n].dx;
                    int ny = y + neighbors[n].dy;
                    int nz = z + neighbors[n].dz;
                    float nval = 1.0f;
                    if (nx >= 0 && nx < gsx && ny >= 0 && ny < gsy && nz >= 0 && nz < gsz) {
                        nval = density_data[nx + gsx * (ny + gsy * nz)];
                    }
                    if (nval <= surf) is_surface = true;
                    gradient += Vector3(neighbors[n].dx, neighbors[n].dy, neighbors[n].dz) * nval;
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
                int zone = zone_data[idx];
                String point_zone = (zone >= 0 && zone < (int)zone_names.size()) ? zone_names[zone] : String();
                if (!zone_matches_fast(point_zone)) continue;

                SampledPoint sampled;
                sampled.normal = normal;
                sampled.zone = zone;
                sampled.zone_name = point_zone;
                sampled.material = (int)material_data[idx];
                sampled.slope = 1.0f - Math::abs(dot_up); // 0 = flat floor, 1 = vertical wall
                sampled.world_pos = Vector3(x, y, z) * voxel_size + normal * (voxel_size * 0.5f);
                local_points.push_back(sampled);
            }
        }
    }
    });

    for (const std::vector<SampledPoint> &points : worker_points) {
        for (const SampledPoint &sampled : points) {
            Dictionary attrs;
            attrs["normal"] = sampled.normal;
            attrs["zone"] = sampled.zone;
            attrs["zone_name"] = sampled.zone_name;
            attrs["material"] = sampled.material;
            attrs["slope"] = sampled.slope;

            Transform3D xform;
            xform.origin = sampled.world_pos;
            point_set->add_raw_point(xform, 1.0f, attrs);
        }
    }

    outputs[0] = context;           // Port 0: Pass-through context
    outputs[1] = point_set;         // Port 1: Generated PointSet
}

} // namespace godot
