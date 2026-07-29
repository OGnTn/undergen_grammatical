#include "undergen_spline_carver_node.h"
#include "density_grid.h"
#include <godot_cpp/classes/resource_loader.hpp>
#include <godot_cpp/core/math.hpp>
#include <godot_cpp/variant/utility_functions.hpp>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

namespace godot {

UnderGenSplineCarverNode::UnderGenSplineCarverNode() {}
UnderGenSplineCarverNode::~UnderGenSplineCarverNode() {}

void UnderGenSplineCarverNode::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_profile_type", "type"), &UnderGenSplineCarverNode::set_profile_type);
    ClassDB::bind_method(D_METHOD("get_profile_type"), &UnderGenSplineCarverNode::get_profile_type);

    ClassDB::bind_method(D_METHOD("set_carve_mode", "mode"), &UnderGenSplineCarverNode::set_carve_mode);
    ClassDB::bind_method(D_METHOD("get_carve_mode"), &UnderGenSplineCarverNode::get_carve_mode);

    ClassDB::bind_method(D_METHOD("set_width", "width"), &UnderGenSplineCarverNode::set_width);
    ClassDB::bind_method(D_METHOD("get_width"), &UnderGenSplineCarverNode::get_width);

    ClassDB::bind_method(D_METHOD("set_height", "height"), &UnderGenSplineCarverNode::set_height);
    ClassDB::bind_method(D_METHOD("get_height"), &UnderGenSplineCarverNode::get_height);

    ClassDB::bind_method(D_METHOD("set_wall_height", "wheight"), &UnderGenSplineCarverNode::set_wall_height);
    ClassDB::bind_method(D_METHOD("get_wall_height"), &UnderGenSplineCarverNode::get_wall_height);

    ClassDB::bind_method(D_METHOD("set_sample_step", "step"), &UnderGenSplineCarverNode::set_sample_step);
    ClassDB::bind_method(D_METHOD("get_sample_step"), &UnderGenSplineCarverNode::get_sample_step);

    ClassDB::bind_method(D_METHOD("set_profile_rotation_deg", "deg"), &UnderGenSplineCarverNode::set_profile_rotation_deg);
    ClassDB::bind_method(D_METHOD("get_profile_rotation_deg"), &UnderGenSplineCarverNode::get_profile_rotation_deg);

    ClassDB::bind_method(D_METHOD("set_custom_profile", "profile"), &UnderGenSplineCarverNode::set_custom_profile);
    ClassDB::bind_method(D_METHOD("get_custom_profile"), &UnderGenSplineCarverNode::get_custom_profile);

    ClassDB::bind_method(D_METHOD("set_curve_resource_path", "path"), &UnderGenSplineCarverNode::set_curve_resource_path);
    ClassDB::bind_method(D_METHOD("get_curve_resource_path"), &UnderGenSplineCarverNode::get_curve_resource_path);

    ClassDB::bind_method(D_METHOD("set_curve", "curve"), &UnderGenSplineCarverNode::set_curve);
    ClassDB::bind_method(D_METHOD("get_curve"), &UnderGenSplineCarverNode::get_curve);

    ClassDB::bind_method(D_METHOD("set_stamp_materials", "stamp"), &UnderGenSplineCarverNode::set_stamp_materials);
    ClassDB::bind_method(D_METHOD("get_stamp_materials"), &UnderGenSplineCarverNode::get_stamp_materials);

    ClassDB::bind_method(D_METHOD("set_floor_material_id", "id"), &UnderGenSplineCarverNode::set_floor_material_id);
    ClassDB::bind_method(D_METHOD("get_floor_material_id"), &UnderGenSplineCarverNode::get_floor_material_id);

    ClassDB::bind_method(D_METHOD("set_wall_material_id", "id"), &UnderGenSplineCarverNode::set_wall_material_id);
    ClassDB::bind_method(D_METHOD("get_wall_material_id"), &UnderGenSplineCarverNode::get_wall_material_id);

    ClassDB::bind_method(D_METHOD("set_ceiling_material_id", "id"), &UnderGenSplineCarverNode::set_ceiling_material_id);
    ClassDB::bind_method(D_METHOD("get_ceiling_material_id"), &UnderGenSplineCarverNode::get_ceiling_material_id);

    ADD_PROPERTY(PropertyInfo(Variant::INT, "profile_type", PROPERTY_HINT_ENUM, "Circle,Horseshoe,Rectangle,Gothic Arch,Custom"), "set_profile_type", "get_profile_type");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "carve_mode", PROPERTY_HINT_ENUM, "Subtract (Carve),Add (Build)"), "set_carve_mode", "get_carve_mode");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "width"), "set_width", "get_width");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "height"), "set_height", "get_height");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "wall_height"), "set_wall_height", "get_wall_height");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "sample_step"), "set_sample_step", "get_sample_step");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "profile_rotation_deg", PROPERTY_HINT_RANGE, "-360,360,1"), "set_profile_rotation_deg", "get_profile_rotation_deg");
    ADD_PROPERTY(PropertyInfo(Variant::PACKED_VECTOR2_ARRAY, "custom_profile"), "set_custom_profile", "get_custom_profile");
    ADD_PROPERTY(PropertyInfo(Variant::STRING, "curve_resource_path"), "set_curve_resource_path", "get_curve_resource_path");
    ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "curve", PROPERTY_HINT_RESOURCE_TYPE, "Curve3D"), "set_curve", "get_curve");

    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "stamp_materials"), "set_stamp_materials", "get_stamp_materials");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "floor_material_id", PROPERTY_HINT_RANGE, "0,255,1"), "set_floor_material_id", "get_floor_material_id");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "wall_material_id", PROPERTY_HINT_RANGE, "0,255,1"), "set_wall_material_id", "get_wall_material_id");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "ceiling_material_id", PROPERTY_HINT_RANGE, "0,255,1"), "set_ceiling_material_id", "get_ceiling_material_id");
}

void UnderGenSplineCarverNode::set_profile_type(ProfileType p_type) { profile_type = p_type; }
UnderGenSplineCarverNode::ProfileType UnderGenSplineCarverNode::get_profile_type() const { return profile_type; }

void UnderGenSplineCarverNode::set_carve_mode(CarveMode p_mode) { carve_mode = p_mode; }
UnderGenSplineCarverNode::CarveMode UnderGenSplineCarverNode::get_carve_mode() const { return carve_mode; }

void UnderGenSplineCarverNode::set_width(float p_width) { width = p_width > 0.1f ? p_width : 0.1f; }
float UnderGenSplineCarverNode::get_width() const { return width; }

void UnderGenSplineCarverNode::set_height(float p_height) { height = p_height > 0.1f ? p_height : 0.1f; }
float UnderGenSplineCarverNode::get_height() const { return height; }

void UnderGenSplineCarverNode::set_wall_height(float p_wheight) { wall_height = p_wheight >= 0.0f ? p_wheight : 0.0f; }
float UnderGenSplineCarverNode::get_wall_height() const { return wall_height; }

void UnderGenSplineCarverNode::set_sample_step(float p_step) { sample_step = p_step > 0.05f ? p_step : 0.05f; }
float UnderGenSplineCarverNode::get_sample_step() const { return sample_step; }

void UnderGenSplineCarverNode::set_profile_rotation_deg(float p_deg) { profile_rotation_deg = p_deg; }
float UnderGenSplineCarverNode::get_profile_rotation_deg() const { return profile_rotation_deg; }

void UnderGenSplineCarverNode::set_custom_profile(const PackedVector2Array &p_profile) { custom_profile = p_profile; }
PackedVector2Array UnderGenSplineCarverNode::get_custom_profile() const { return custom_profile; }

void UnderGenSplineCarverNode::set_curve_resource_path(const String &p_path) { curve_resource_path = p_path; }
String UnderGenSplineCarverNode::get_curve_resource_path() const { return curve_resource_path; }

void UnderGenSplineCarverNode::set_curve(const Ref<Curve3D> &p_curve) { curve = p_curve; }
Ref<Curve3D> UnderGenSplineCarverNode::get_curve() const { return curve; }

void UnderGenSplineCarverNode::set_stamp_materials(bool p_stamp) { stamp_materials = p_stamp; }
bool UnderGenSplineCarverNode::get_stamp_materials() const { return stamp_materials; }

void UnderGenSplineCarverNode::set_floor_material_id(int p_id) { floor_material_id = Math::clamp(p_id, 0, 255); }
int UnderGenSplineCarverNode::get_floor_material_id() const { return floor_material_id; }

void UnderGenSplineCarverNode::set_wall_material_id(int p_id) { wall_material_id = Math::clamp(p_id, 0, 255); }
int UnderGenSplineCarverNode::get_wall_material_id() const { return wall_material_id; }

void UnderGenSplineCarverNode::set_ceiling_material_id(int p_id) { ceiling_material_id = Math::clamp(p_id, 0, 255); }
int UnderGenSplineCarverNode::get_ceiling_material_id() const { return ceiling_material_id; }

float UnderGenSplineCarverNode::_evaluate_2d_sdf(float u, float v) const {
    float half_w = width * 0.5f;

    switch (profile_type) {
        case PROFILE_CIRCLE: {
            float r = half_w;
            return std::sqrt(u * u + v * v) - r;
        }
        case PROFILE_RECTANGLE: {
            float dx = std::abs(u) - half_w;
            float dy = std::abs(v - height * 0.5f) - height * 0.5f;
            return std::max(dx, dy);
        }
        case PROFILE_HORSESHOE: {
            float wheight = (wall_height < height) ? wall_height : height * 0.6f;
            if (v <= wheight) {
                float dx = std::abs(u) - half_w;
                float dy = (v < 0.0f) ? -v : (v - wheight);
                return std::max(dx, dy);
            } else {
                float dv = v - wheight;
                float dist_arch = std::sqrt(u * u + dv * dv) - half_w;
                float dx = std::abs(u) - half_w;
                return std::max(dist_arch, dx);
            }
        }
        case PROFILE_GOTHIC_ARCH: {
            float wheight = (wall_height < height) ? wall_height : height * 0.5f;
            if (v <= wheight) {
                float dx = std::abs(u) - half_w;
                float dy = (v < 0.0f) ? -v : (v - wheight);
                return std::max(dx, dy);
            } else {
                float dv = v - wheight;
                float dist_left = std::sqrt((u + half_w) * (u + half_w) + dv * dv) - width;
                float dist_right = std::sqrt((u - half_w) * (u - half_w) + dv * dv) - width;
                return std::max(dist_left, dist_right);
            }
        }
        case PROFILE_CUSTOM: {
            if (custom_profile.size() < 3) {
                return std::sqrt(u * u + v * v) - half_w;
            }
            Vector2 pt(u, v);
            float min_dist_sq = 1e9f;
            bool inside = false;
            int n = custom_profile.size();
            for (int i = 0, j = n - 1; i < n; j = i++) {
                Vector2 p1 = custom_profile[i];
                Vector2 p2 = custom_profile[j];

                Vector2 seg = p2 - p1;
                float l2 = seg.length_squared();
                float t = (l2 > 0.0f) ? std::max(0.0f, std::min(1.0f, (pt - p1).dot(seg) / l2)) : 0.0f;
                Vector2 proj = p1 + t * seg;
                float dist_sq = (pt - proj).length_squared();
                if (dist_sq < min_dist_sq) min_dist_sq = dist_sq;

                if (((p1.y > pt.y) != (p2.y > pt.y)) &&
                    (pt.x < (p2.x - p1.x) * (pt.y - p1.y) / (p2.y - p1.y + 1e-6f) + p1.x)) {
                    inside = !inside;
                }
            }
            float dist = std::sqrt(min_dist_sq);
            return inside ? -dist : dist;
        }
    }
    return std::sqrt(u * u + v * v) - half_w;
}

struct CurveFrame {
    Vector3 pos;
    Vector3 tangent;
    Vector3 normal;
    Vector3 binormal;
};

void UnderGenSplineCarverNode::_carve_single_curve(Ref<DensityGrid> grid, const Ref<Curve3D> &active_curve) {
    if (active_curve.is_null() || active_curve->get_point_count() < 2) return;

    float baked_len = active_curve->get_baked_length();
    if (baked_len <= 0.1f) return;

    int sample_count = std::max(2, (int)std::ceil(baked_len / sample_step));
    std::vector<CurveFrame> frames(sample_count);

    Vector3 p0 = active_curve->sample_baked(0.0f);
    Vector3 p1 = active_curve->sample_baked(Math::min(0.1f, baked_len));
    Vector3 t0 = (p1 - p0).normalized();
    if (t0.length_squared() < 0.001f) t0 = Vector3(0, 0, 1);

    Vector3 up(0, 1, 0);
    if (std::abs(t0.dot(up)) > 0.95f) up = Vector3(1, 0, 0);

    // Initial frame: b0 = UP, n0 = RIGHT
    Vector3 b0 = (up - t0 * t0.dot(up)).normalized();
    Vector3 n0 = t0.cross(b0).normalized();

    frames[0] = { p0, t0, n0, b0 };

    for (int i = 1; i < sample_count; ++i) {
        float distance = (float)i / (sample_count - 1) * baked_len;
        Vector3 pos = active_curve->sample_baked(distance);
        Vector3 prev_pos = frames[i - 1].pos;
        Vector3 tan = (pos - prev_pos).normalized();
        if (tan.length_squared() < 0.001f) tan = frames[i - 1].tangent;

        Vector3 v1 = pos - prev_pos;
        float c1 = v1.dot(v1);
        Vector3 r_n = frames[i - 1].normal - (v1 * (2.0f * v1.dot(frames[i - 1].normal) / c1));
        Vector3 r_t = frames[i - 1].tangent - (v1 * (2.0f * v1.dot(frames[i - 1].tangent) / c1));

        Vector3 v2 = tan - r_t;
        float c2 = v2.dot(v2);
        Vector3 norm = (c2 > 0.0001f) ? (r_n - (v2 * (2.0f * v2.dot(r_n) / c2))).normalized() : r_n.normalized();
        Vector3 binorm = tan.cross(norm).normalized();

        frames[i] = { pos, tan, norm, binorm };
    }

    Vector3 min_bounds(1e9f, 1e9f, 1e9f);
    Vector3 max_bounds(-1e9f, -1e9f, -1e9f);
    float margin = std::max(width, height) * 0.75f + 2.0f;

    for (const auto &f : frames) {
        min_bounds.x = Math::min(min_bounds.x, f.pos.x - margin);
        min_bounds.y = Math::min(min_bounds.y, f.pos.y - margin);
        min_bounds.z = Math::min(min_bounds.z, f.pos.z - margin);

        max_bounds.x = Math::max(max_bounds.x, f.pos.x + margin);
        max_bounds.y = Math::max(max_bounds.y, f.pos.y + margin);
        max_bounds.z = Math::max(max_bounds.z, f.pos.z + margin);
    }

    Vector3i grid_dims = grid->get_grid_dimensions();
    Vector3i min_v(
        Math::clamp((int)Math::floor(min_bounds.x), 0, grid_dims.x - 1),
        Math::clamp((int)Math::floor(min_bounds.y), 0, grid_dims.y - 1),
        Math::clamp((int)Math::floor(min_bounds.z), 0, grid_dims.z - 1)
    );
    Vector3i max_v(
        Math::clamp((int)Math::ceil(max_bounds.x), 0, grid_dims.x - 1),
        Math::clamp((int)Math::ceil(max_bounds.y), 0, grid_dims.y - 1),
        Math::clamp((int)Math::ceil(max_bounds.z), 0, grid_dims.z - 1)
    );

    float rot_rad = Math::deg_to_rad(profile_rotation_deg);
    float cos_r = std::cos(rot_rad);
    float sin_r = std::sin(rot_rad);

    for (int z = min_v.z; z <= max_v.z; ++z) {
        for (int y = min_v.y; y <= max_v.y; ++y) {
            for (int x = min_v.x; x <= max_v.x; ++x) {
                Vector3 voxel_pos((float)x, (float)y, (float)z);

                int closest_idx = 0;
                float min_dist_sq = 1e9f;
                for (size_t i = 0; i < frames.size(); ++i) {
                    float d_sq = (voxel_pos - frames[i].pos).length_squared();
                    if (d_sq < min_dist_sq) {
                        min_dist_sq = d_sq;
                        closest_idx = (int)i;
                    }
                }

                const CurveFrame &frame = frames[closest_idx];
                Vector3 offset = voxel_pos - frame.pos;

                float u_raw = offset.dot(frame.normal);
                float v_raw = offset.dot(frame.binormal);

                float u = u_raw * cos_r - v_raw * sin_r;
                float v = u_raw * sin_r + v_raw * cos_r;

                float sdf = _evaluate_2d_sdf(u, v);
                Vector3i cell_pos(x, y, z);
                float current_density = grid->get_cell(cell_pos);

                if (carve_mode == CARVE_SUBTRACT) {
                    if (sdf < 0.0f) {
                        grid->set_cell(cell_pos, -1.0f);
                    } else if (sdf < 1.0f) {
                        float target = Math::lerp(-1.0f, current_density, sdf);
                        grid->set_cell(cell_pos, Math::min(current_density, target));
                    }

                    if (stamp_materials && sdf <= 0.5f) {
                        if (v <= 0.5f && floor_material_id > 0) {
                            grid->set_material_id(cell_pos, floor_material_id);
                        } else if (v >= height - 0.8f && ceiling_material_id > 0) {
                            grid->set_material_id(cell_pos, ceiling_material_id);
                        } else if (wall_material_id > 0) {
                            grid->set_material_id(cell_pos, wall_material_id);
                        }
                    }
                } else if (carve_mode == CARVE_ADD) {
                    if (sdf < 0.0f) {
                        grid->set_cell(cell_pos, 1.0f);
                        if (stamp_materials && wall_material_id > 0) {
                            grid->set_material_id(cell_pos, wall_material_id);
                        }
                    }
                }
            }
        }
    }
}

void UnderGenSplineCarverNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    Dictionary context = inputs.get(0, Dictionary());
    if (context.is_empty()) return;

    Ref<DensityGrid> grid = context.get("grid", Ref<DensityGrid>());
    if (grid.is_null()) {
        outputs[0] = context;
        return;
    }

    // 1. Check for explicitly assigned single Curve3D
    Ref<Curve3D> active_curve = curve;
    if (active_curve.is_null() && !curve_resource_path.is_empty()) {
        active_curve = ResourceLoader::get_singleton()->load(curve_resource_path);
    }
    if (active_curve.is_null() && context.has("curve")) {
        active_curve = context["curve"];
    }

    if (active_curve.is_valid() && active_curve->get_point_count() >= 2) {
        _carve_single_curve(grid, active_curve);
        outputs[0] = context;
        return;
    }

    // 2. Otherwise, fallback to generating curves between rooms in context["rooms"] and context["edges"]
    Array rooms_arr = context.get("rooms", Array());
    Array edges_arr = context.get("edges", Array());

    if (rooms_arr.size() >= 2 && !edges_arr.is_empty()) {
        std::map<String, Vector3> room_positions;
        for (int i = 0; i < rooms_arr.size(); ++i) {
            Dictionary r = rooms_arr[i];
            String id = r.get("id", "");
            Vector3i pos = r.get("position", Vector3i());
            Vector3i size = r.get("size", Vector3i());
            Vector3 center = Vector3(pos) + Vector3(size) * 0.5f;
            room_positions[id] = center;
        }

        for (int i = 0; i < edges_arr.size(); ++i) {
            Dictionary edge = edges_arr[i];
            String from_id = edge.get("from", "");
            String to_id = edge.get("to", "");

            if (room_positions.count(from_id) && room_positions.count(to_id)) {
                Vector3 p_from = room_positions[from_id];
                Vector3 p_to = room_positions[to_id];

                Ref<Curve3D> edge_curve;
                edge_curve.instantiate();
                edge_curve->add_point(p_from);

                Vector3 mid_point = Vector3(p_to.x, (p_from.y + p_to.y) * 0.5f, p_from.z);
                edge_curve->add_point(mid_point);
                edge_curve->add_point(p_to);

                _carve_single_curve(grid, edge_curve);
            }
        }
    }

    outputs[0] = context;
}

} // namespace godot
