#include "undergen_outdoor_placer_node.h"
#include "level_gen_data.h"
#include <godot_cpp/variant/utility_functions.hpp>
#include <godot_cpp/core/math.hpp>
#include <godot_cpp/classes/file_access.hpp>
#include "ogt_vox.h"
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>

namespace godot {

static float get_dist_to_box(float x, float z, const Vector3i &pos, const Vector3i &size) {
    float min_x = pos.x;
    float max_x = pos.x + size.x;
    float min_z = pos.z;
    float max_z = pos.z + size.z;
    float dx = std::max({min_x - x, 0.0f, x - max_x});
    float dz = std::max({min_z - z, 0.0f, z - max_z});
    return std::sqrt(dx * dx + dz * dz);
}

static float get_dist_to_segment(float px, float pz, float ax, float az, float bx, float bz, float &out_t) {
    float abx = bx - ax;
    float abz = bz - az;
    float apx = px - ax;
    float apz = pz - az;
    float ab_len_sq = abx * abx + abz * abz;
    if (ab_len_sq < 0.0001f) {
        out_t = 0.0f;
        return std::sqrt(apx * apx + apz * apz);
    }
    float t = (apx * abx + apz * abz) / ab_len_sq;
    t = Math::clamp(t, 0.0f, 1.0f);
    out_t = t;
    float projx = ax + t * abx;
    float projz = az + t * abz;
    float dx = px - projx;
    float dz = pz - projz;
    return std::sqrt(dx * dx + dz * dz);
}

UnderGenOutdoorPlacerNode::UnderGenOutdoorPlacerNode() {
    rng.instantiate();
    noise_generator.instantiate();
}

UnderGenOutdoorPlacerNode::~UnderGenOutdoorPlacerNode() {}

void UnderGenOutdoorPlacerNode::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_grid_size_x", "grid_size_x"), &UnderGenOutdoorPlacerNode::set_grid_size_x);
    ClassDB::bind_method(D_METHOD("get_grid_size_x"), &UnderGenOutdoorPlacerNode::get_grid_size_x);
    ClassDB::bind_method(D_METHOD("set_grid_size_y", "grid_size_y"), &UnderGenOutdoorPlacerNode::set_grid_size_y);
    ClassDB::bind_method(D_METHOD("get_grid_size_y"), &UnderGenOutdoorPlacerNode::get_grid_size_y);
    ClassDB::bind_method(D_METHOD("set_grid_size_z", "grid_size_z"), &UnderGenOutdoorPlacerNode::set_grid_size_z);
    ClassDB::bind_method(D_METHOD("get_grid_size_z"), &UnderGenOutdoorPlacerNode::get_grid_size_z);
    ClassDB::bind_method(D_METHOD("set_surface_threshold", "threshold"), &UnderGenOutdoorPlacerNode::set_surface_threshold);
    ClassDB::bind_method(D_METHOD("get_surface_threshold"), &UnderGenOutdoorPlacerNode::get_surface_threshold);

    ClassDB::bind_method(D_METHOD("set_margin_x", "margin_x"), &UnderGenOutdoorPlacerNode::set_margin_x);
    ClassDB::bind_method(D_METHOD("get_margin_x"), &UnderGenOutdoorPlacerNode::get_margin_x);
    ClassDB::bind_method(D_METHOD("set_margin_y", "margin_y"), &UnderGenOutdoorPlacerNode::set_margin_y);
    ClassDB::bind_method(D_METHOD("get_margin_y"), &UnderGenOutdoorPlacerNode::get_margin_y);
    ClassDB::bind_method(D_METHOD("set_margin_z", "margin_z"), &UnderGenOutdoorPlacerNode::set_margin_z);
    ClassDB::bind_method(D_METHOD("get_margin_z"), &UnderGenOutdoorPlacerNode::get_margin_z);

    ClassDB::bind_method(D_METHOD("set_spread_ratio", "ratio"), &UnderGenOutdoorPlacerNode::set_spread_ratio);
    ClassDB::bind_method(D_METHOD("get_spread_ratio"), &UnderGenOutdoorPlacerNode::get_spread_ratio);

    ClassDB::bind_method(D_METHOD("set_base_height", "height"), &UnderGenOutdoorPlacerNode::set_base_height);
    ClassDB::bind_method(D_METHOD("get_base_height"), &UnderGenOutdoorPlacerNode::get_base_height);
    ClassDB::bind_method(D_METHOD("set_mountain_height", "height"), &UnderGenOutdoorPlacerNode::set_mountain_height);
    ClassDB::bind_method(D_METHOD("get_mountain_height"), &UnderGenOutdoorPlacerNode::get_mountain_height);
    ClassDB::bind_method(D_METHOD("set_slope_width", "width"), &UnderGenOutdoorPlacerNode::set_slope_width);
    ClassDB::bind_method(D_METHOD("get_slope_width"), &UnderGenOutdoorPlacerNode::get_slope_width);
    ClassDB::bind_method(D_METHOD("set_path_width", "width"), &UnderGenOutdoorPlacerNode::set_path_width);
    ClassDB::bind_method(D_METHOD("get_path_width"), &UnderGenOutdoorPlacerNode::get_path_width);
    ClassDB::bind_method(D_METHOD("set_noise_intensity", "intensity"), &UnderGenOutdoorPlacerNode::set_noise_intensity);
    ClassDB::bind_method(D_METHOD("get_noise_intensity"), &UnderGenOutdoorPlacerNode::get_noise_intensity);
    ClassDB::bind_method(D_METHOD("set_noise_generator", "noise"), &UnderGenOutdoorPlacerNode::set_noise_generator);
    ClassDB::bind_method(D_METHOD("get_noise_generator"), &UnderGenOutdoorPlacerNode::get_noise_generator);

    ADD_PROPERTY(PropertyInfo(Variant::INT, "grid_size_x"), "set_grid_size_x", "get_grid_size_x");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "grid_size_y"), "set_grid_size_y", "get_grid_size_y");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "grid_size_z"), "set_grid_size_z", "get_grid_size_z");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "margin_x"), "set_margin_x", "get_margin_x");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "margin_y"), "set_margin_y", "get_margin_y");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "margin_z"), "set_margin_z", "get_margin_z");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "surface_threshold"), "set_surface_threshold", "get_surface_threshold");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "spread_ratio", PROPERTY_HINT_RANGE, "0.0,1.0,0.01"), "set_spread_ratio", "get_spread_ratio");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "base_height"), "set_base_height", "get_base_height");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "mountain_height"), "set_mountain_height", "get_mountain_height");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "slope_width"), "set_slope_width", "get_slope_width");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "path_width"), "set_path_width", "get_path_width");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "noise_intensity"), "set_noise_intensity", "get_noise_intensity");
    ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "noise_generator", PROPERTY_HINT_RESOURCE_TYPE, "FastNoiseLite"), "set_noise_generator", "get_noise_generator");
}

void UnderGenOutdoorPlacerNode::set_grid_size_x(int p_x) { grid_size_x = p_x; }
int UnderGenOutdoorPlacerNode::get_grid_size_x() const { return grid_size_x; }
void UnderGenOutdoorPlacerNode::set_grid_size_y(int p_y) { grid_size_y = p_y; }
int UnderGenOutdoorPlacerNode::get_grid_size_y() const { return grid_size_y; }
void UnderGenOutdoorPlacerNode::set_grid_size_z(int p_z) { grid_size_z = p_z; }
int UnderGenOutdoorPlacerNode::get_grid_size_z() const { return grid_size_z; }
void UnderGenOutdoorPlacerNode::set_surface_threshold(float p_threshold) { surface_threshold = p_threshold; }
float UnderGenOutdoorPlacerNode::get_surface_threshold() const { return surface_threshold; }

void UnderGenOutdoorPlacerNode::set_margin_x(int p_x) { margin_x = p_x; }
int UnderGenOutdoorPlacerNode::get_margin_x() const { return margin_x; }
void UnderGenOutdoorPlacerNode::set_margin_y(int p_y) { margin_y = p_y; }
int UnderGenOutdoorPlacerNode::get_margin_y() const { return margin_y; }
void UnderGenOutdoorPlacerNode::set_margin_z(int p_z) { margin_z = p_z; }
int UnderGenOutdoorPlacerNode::get_margin_z() const { return margin_z; }

void UnderGenOutdoorPlacerNode::set_spread_ratio(float p_ratio) { spread_ratio = p_ratio; }
float UnderGenOutdoorPlacerNode::get_spread_ratio() const { return spread_ratio; }

void UnderGenOutdoorPlacerNode::set_base_height(int p_height) { base_height = p_height; }
int UnderGenOutdoorPlacerNode::get_base_height() const { return base_height; }
void UnderGenOutdoorPlacerNode::set_mountain_height(int p_height) { mountain_height = p_height; }
int UnderGenOutdoorPlacerNode::get_mountain_height() const { return mountain_height; }
void UnderGenOutdoorPlacerNode::set_slope_width(float p_width) { slope_width = p_width; }
float UnderGenOutdoorPlacerNode::get_slope_width() const { return slope_width; }
void UnderGenOutdoorPlacerNode::set_path_width(float p_width) { path_width = p_width; }
float UnderGenOutdoorPlacerNode::get_path_width() const { return path_width; }
void UnderGenOutdoorPlacerNode::set_noise_intensity(float p_intensity) { noise_intensity = p_intensity; }
float UnderGenOutdoorPlacerNode::get_noise_intensity() const { return noise_intensity; }
void UnderGenOutdoorPlacerNode::set_noise_generator(const Ref<FastNoiseLite> &p_noise) { noise_generator = p_noise; }
Ref<FastNoiseLite> UnderGenOutdoorPlacerNode::get_noise_generator() const { return noise_generator; }

void UnderGenOutdoorPlacerNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    int64_t seed = inputs.get(0, 12345);
    Dictionary logical_graph = inputs.get(1, Dictionary());
    
    rng->set_seed(seed);

    Ref<DensityGrid> grid;
    grid.instantiate();
    // Default open/empty (air = -1.0f) for outdoor level
    grid->initialize_grid(grid_size_x, grid_size_y, grid_size_z, -1.0f);
    grid->set_surface_threshold(surface_threshold);

    Array nodes = logical_graph.get("nodes", Array());
    Array edges = logical_graph.get("edges", Array());

    std::vector<ResolvedRoom> processing_rooms;
    for (int i = 0; i < nodes.size(); ++i) {
        Dictionary node = nodes[i];
        ResolvedRoom r;
        r.id = node.get("id", "room");
        r.type = node.get("type", "generic");
        r.vox_path = node.get("vox_path", "");
        r.exclude_from_smoothing = node.get("exclude_from_smoothing", false);
        r.exclude_from_warping = node.get("exclude_from_warping", false);
        
        Vector3i min_s = node.get("min_size", Vector3i(5, 5, 5));
        Vector3i max_s = node.get("max_size", Vector3i(10, 10, 10));
        bool size_loaded = false;

        if (!r.vox_path.is_empty()) {
            Ref<FileAccess> file = FileAccess::open(r.vox_path, FileAccess::READ);
            if (file.is_valid()) {
                uint64_t len = file->get_length();
                PackedByteArray buf = file->get_buffer(len);
                const ogt_vox_scene* scene = ogt_vox_read_scene(buf.ptr(), (uint32_t)len);
                if (scene) {
                    int scene_min_x = INT_MAX, scene_max_x = INT_MIN;
                    int scene_min_y = INT_MAX, scene_max_y = INT_MIN;
                    int scene_min_z = INT_MAX, scene_max_z = INT_MIN;
                    bool has_visible_instances = false;

                    for (uint32_t j = 0; j < scene->num_instances; ++j) {
                        const ogt_vox_instance& inst = scene->instances[j];
                        if (inst.hidden) continue;
                        const ogt_vox_model* model = scene->models[inst.model_index];
                        if (!model) continue;

                        has_visible_instances = true;

                        int offset_x = (int)inst.transform.m30;
                        int offset_y = (int)inst.transform.m31;
                        int offset_z = (int)inst.transform.m32;

                        int inst_min_x = offset_x;
                        int inst_max_x = offset_x + model->size_x;
                        int inst_min_y = offset_z; // MV Z -> Godot Y
                        int inst_max_y = offset_z + model->size_z;
                        int inst_min_z = offset_y; // MV Y -> Godot Z
                        int inst_max_z = offset_y + model->size_y;

                        if (inst_min_x < scene_min_x) scene_min_x = inst_min_x;
                        if (inst_max_x > scene_max_x) scene_max_x = inst_max_x;
                        if (inst_min_y < scene_min_y) scene_min_y = inst_min_y;
                        if (inst_max_y > scene_max_y) scene_max_y = inst_max_y;
                        if (inst_min_z < scene_min_z) scene_min_z = inst_min_z;
                        if (inst_max_z > scene_max_z) scene_max_z = inst_max_z;
                    }

                    if (has_visible_instances) {
                        r.size = Vector3i(
                            scene_max_x - scene_min_x,
                            scene_max_y - scene_min_y,
                            scene_max_z - scene_min_z
                        );
                        size_loaded = true;
                    }
                    ogt_vox_destroy_scene(scene);
                }
            }
        }

        if (!size_loaded) {
            r.size = Vector3i(
                rng->randi_range(min_s.x, max_s.x),
                rng->randi_range(min_s.y, max_s.y),
                rng->randi_range(min_s.z, max_s.z)
            );
        }
        
        Dictionary constraints = node.get("constraints", Dictionary());
        Vector3i grid_center = grid->get_grid_dimensions() / 2;
        if (constraints.has("fixed_pos")) {
             Vector3 fix_v = constraints["fixed_pos"];
             r.is_fixed = true;
             r.fixed_position = Vector3(grid_center) + fix_v;
             r.position = Vector3i(r.fixed_position) - (r.size / 2);
        } else {
             Vector3i max_pos = grid->get_grid_dimensions() - r.size;
             r.position = Vector3i(
                 rng->randi_range(0, MAX(0, max_pos.x)),
                 base_height,
                 rng->randi_range(0, MAX(0, max_pos.z))
             );
        }
        
        processing_rooms.push_back(r);
    }

    // Run Collision-Free Random Placement
    if (!processing_rooms.empty()) {
        int safe_margin_x = Math::clamp(margin_x, 1, Math::max(1, (grid_size_x - 6) / 2));
        int safe_margin_z = Math::clamp(margin_z, 1, Math::max(1, (grid_size_z - 6) / 2));

        Vector3i grid_center = grid->get_grid_dimensions() / 2;

        std::vector<size_t> fixed_room_indices;
        std::vector<size_t> unfixed_room_indices;
        for (size_t i = 0; i < processing_rooms.size(); ++i) {
            if (processing_rooms[i].is_fixed) {
                fixed_room_indices.push_back(i);
            } else {
                unfixed_room_indices.push_back(i);
            }
        }

        std::sort(unfixed_room_indices.begin(), unfixed_room_indices.end(), [&](size_t a, size_t b) {
            Vector3i sa = processing_rooms[a].size;
            Vector3i sb = processing_rooms[b].size;
            return (sa.x * sa.y * sa.z) > (sb.x * sb.y * sb.z);
        });

        std::vector<size_t> placed_room_indices = fixed_room_indices;

        for (size_t idx : unfixed_room_indices) {
            ResolvedRoom &room = processing_rooms[idx];

            Vector3i min_bounds(safe_margin_x, base_height, safe_margin_z);
            Vector3i max_bounds(
                grid_size_x - safe_margin_x - room.size.x,
                base_height,
                grid_size_z - safe_margin_z - room.size.z
            );

            std::vector<Vector3i> candidates;
            int current_spacing = 8;

            while (candidates.empty() && current_spacing >= 0) {
                for (int attempt = 0; attempt < 1000; ++attempt) {
                    int rx = (max_bounds.x > min_bounds.x) ? rng->randi_range(min_bounds.x, max_bounds.x) : min_bounds.x;
                    int ry = base_height; 
                    int rz = (max_bounds.z > min_bounds.z) ? rng->randi_range(min_bounds.z, max_bounds.z) : min_bounds.z;
                    Vector3i proposed_pos(rx, ry, rz);

                    bool overlaps = false;
                    for (size_t p_idx : placed_room_indices) {
                        const ResolvedRoom &placed = processing_rooms[p_idx];
                        bool overlap_x = (proposed_pos.x - current_spacing < placed.position.x + placed.size.x) &&
                                         (proposed_pos.x + room.size.x + current_spacing > placed.position.x);
                        bool overlap_y = (proposed_pos.y - current_spacing < placed.position.y + placed.size.y) &&
                                         (proposed_pos.y + room.size.y + current_spacing > placed.position.y);
                        bool overlap_z = (proposed_pos.z - current_spacing < placed.position.z + placed.size.z) &&
                                         (proposed_pos.z + room.size.z + current_spacing > placed.position.z);
                        if (overlap_x && overlap_y && overlap_z) {
                            overlaps = true;
                            break;
                        }
                    }

                    if (!overlaps) {
                        candidates.push_back(proposed_pos);
                        if (candidates.size() >= 30) {
                            break;
                        }
                    }
                }

                if (candidates.empty()) {
                    current_spacing -= 2;
                }
            }

            Vector3i chosen_pos;
            if (!candidates.empty()) {
                std::sort(candidates.begin(), candidates.end(), [&](const Vector3i &a, const Vector3i &b) {
                    Vector3 center_a(a.x + room.size.x / 2.0f, a.y + room.size.y / 2.0f, a.z + room.size.z / 2.0f);
                    Vector3 center_b(b.x + room.size.x / 2.0f, b.y + room.size.y / 2.0f, b.z + room.size.z / 2.0f);
                    Vector3 center_grid(grid_center.x, grid_center.y, grid_center.z);
                    return center_a.distance_to(center_grid) < center_b.distance_to(center_grid);
                });

                float t = Math::abs(spread_ratio - 0.5f) * 2.0f;
                int chosen_idx = 0;
                if (spread_ratio < 0.5f) {
                    int rand_idx = rng->randi() % candidates.size();
                    chosen_idx = (int)Math::round(Math::lerp((float)rand_idx, 0.0f, t));
                } else {
                    int rand_idx = rng->randi() % candidates.size();
                    chosen_idx = (int)Math::round(Math::lerp((float)rand_idx, (float)(candidates.size() - 1), t));
                }
                chosen_pos = candidates[chosen_idx];
            } else {
                int rx = (max_bounds.x > min_bounds.x) ? rng->randi_range(min_bounds.x, max_bounds.x) : min_bounds.x;
                int ry = base_height;
                int rz = (max_bounds.z > min_bounds.z) ? rng->randi_range(min_bounds.z, max_bounds.z) : min_bounds.z;
                chosen_pos = Vector3i(rx, ry, rz);
            }

            room.position = chosen_pos;
            placed_room_indices.push_back(idx);
        }
    }

    std::map<String, int> id_to_index;
    for (size_t i = 0; i < processing_rooms.size(); ++i) {
        id_to_index[processing_rooms[i].id] = i;
    }

    struct PathSegment {
        Vector2 start_pos;
        Vector2 end_pos;
        float start_height;
        float end_height;
    };

    std::vector<PathSegment> segments;
    for (int i = 0; i < edges.size(); ++i) {
        Dictionary e_dict = edges[i];
        String from = e_dict.get("from", "");
        String to = e_dict.get("to", "");
        if (id_to_index.count(from) && id_to_index.count(to)) {
            const ResolvedRoom &r_from = processing_rooms[id_to_index[from]];
            const ResolvedRoom &r_to = processing_rooms[id_to_index[to]];
            PathSegment seg;
            seg.start_pos = Vector2(r_from.center().x, r_from.center().z);
            seg.end_pos = Vector2(r_to.center().x, r_to.center().z);
            seg.start_height = (float)r_from.position.y;
            seg.end_height = (float)r_to.position.y;
            segments.push_back(seg);
        }
    }

    float slope_w = Math::max(0.1f, slope_width);
    float path_w_half = path_width * 0.5f;

    for (int z = 0; z < grid_size_z; ++z) {
        for (int x = 0; x < grid_size_x; ++x) {
            float mountain_h = (float)mountain_height;
            if (noise_generator.is_valid()) {
                mountain_h += noise_generator->get_noise_2d((float)x, (float)z) * noise_intensity;
            }
            mountain_h = Math::clamp(mountain_h, 0.0f, (float)grid_size_y - 1.0f);

            float min_d_room = INFINITY;
            float target_h_room = (float)base_height;
            for (const auto &room : processing_rooms) {
                float d = get_dist_to_box((float)x, (float)z, room.position, room.size);
                if (d < min_d_room) {
                    min_d_room = d;
                    target_h_room = (float)room.position.y;
                }
            }

            float min_d_path = INFINITY;
            float target_h_path = (float)base_height;
            for (const auto &seg : segments) {
                float t = 0.0f;
                float d = get_dist_to_segment((float)x, (float)z, seg.start_pos.x, seg.start_pos.y, seg.end_pos.x, seg.end_pos.y, t);
                float d_adj = Math::max(0.0f, d - path_w_half);
                if (d_adj < min_d_path) {
                    min_d_path = d_adj;
                    target_h_path = Math::lerp(seg.start_height, seg.end_height, t);
                }
            }

            float d_feature = Math::min(min_d_room, min_d_path);
            float target_h_feature = (min_d_room < min_d_path) ? target_h_room : target_h_path;

            float final_h = mountain_h;
            if (d_feature <= 0.0f) {
                final_h = target_h_feature;
            } else if (d_feature < slope_w) {
                float t = d_feature / slope_w;
                final_h = Math::lerp(target_h_feature, mountain_h, t);
            }

            for (int y = 0; y < grid_size_y; ++y) {
                Vector3i cell_pos(x, y, z);
                if ((float)y < final_h) {
                    grid->set_cell(cell_pos, 1.0f); // Solid ground
                } else {
                    grid->set_cell(cell_pos, -1.0f); // Empty air
                }
            }
        }
    }

    // Set Zone IDs for room volumes
    for (size_t i = 0; i < processing_rooms.size(); ++i) {
        const auto &room = processing_rooms[i];
        int z_id = grid->register_zone_name(room.type);
        if (z_id > 0) {
            for (int rx = room.position.x; rx < room.position.x + room.size.x; ++rx) {
                for (int ry = room.position.y; ry < room.position.y + room.size.y; ++ry) {
                    for (int rz = room.position.z; rz < room.position.z + room.size.z; ++rz) {
                        Vector3i pos(rx, ry, rz);
                        if (grid->is_valid_position(pos)) {
                            grid->set_zone_at(pos, z_id);
                        }
                    }
                }
            }
        }
    }

    Array placed_rooms_array;
    for (const auto &room : processing_rooms) {
        Dictionary r_dict;
        r_dict["id"] = room.id;
        r_dict["type"] = room.type;
        r_dict["position"] = room.position;
        r_dict["size"] = room.size;
        r_dict["vox_path"] = room.vox_path;
        r_dict["exclude_from_smoothing"] = room.exclude_from_smoothing;
        r_dict["exclude_from_warping"] = room.exclude_from_warping;
        r_dict["center"] = room.center();
        placed_rooms_array.append(r_dict);
    }

    Dictionary context;
    context["grid"] = grid;
    context["rooms"] = placed_rooms_array;
    context["edges"] = edges;
    context["seed"] = seed;

    outputs[0] = context;
}

Dictionary UnderGenOutdoorPlacerNode::get_pipeline_input_defaults(const Dictionary &global_inputs) const {
    Dictionary defaults;
    if (global_inputs.has(0)) defaults[0] = global_inputs[0];
    return defaults;
}

} // namespace godot
