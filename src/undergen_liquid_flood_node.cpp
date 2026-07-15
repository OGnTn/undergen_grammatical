#include "undergen_liquid_flood_node.h"
#include "density_grid.h"
#include "undergen_point_set.h"
#include <godot_cpp/core/math.hpp>
#include <godot_cpp/variant/utility_functions.hpp>
#include <map>
#include <vector>
#include <queue>
#include <algorithm>

namespace godot {

UnderGenLiquidFloodNode::UnderGenLiquidFloodNode() {
    rng.instantiate();
    rng->randomize();
}

UnderGenLiquidFloodNode::~UnderGenLiquidFloodNode() {}

void UnderGenLiquidFloodNode::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_liquid_material_id", "id"), &UnderGenLiquidFloodNode::set_liquid_material_id);
    ClassDB::bind_method(D_METHOD("get_liquid_material_id"), &UnderGenLiquidFloodNode::get_liquid_material_id);
    ClassDB::bind_method(D_METHOD("set_zone_flood_levels", "levels"), &UnderGenLiquidFloodNode::set_zone_flood_levels);
    ClassDB::bind_method(D_METHOD("get_zone_flood_levels"), &UnderGenLiquidFloodNode::get_zone_flood_levels);
    ClassDB::bind_method(D_METHOD("set_pool_radius", "radius"), &UnderGenLiquidFloodNode::set_pool_radius);
    ClassDB::bind_method(D_METHOD("get_pool_radius"), &UnderGenLiquidFloodNode::get_pool_radius);
    ClassDB::bind_method(D_METHOD("set_pool_depth", "depth"), &UnderGenLiquidFloodNode::set_pool_depth);
    ClassDB::bind_method(D_METHOD("get_pool_depth"), &UnderGenLiquidFloodNode::get_pool_depth);
    ClassDB::bind_method(D_METHOD("set_voxel_size", "size"), &UnderGenLiquidFloodNode::set_voxel_size);
    ClassDB::bind_method(D_METHOD("get_voxel_size"), &UnderGenLiquidFloodNode::get_voxel_size);

    // Bind Flow & Basin properties
    ClassDB::bind_method(D_METHOD("set_flow_mode", "mode"), &UnderGenLiquidFloodNode::set_flow_mode);
    ClassDB::bind_method(D_METHOD("get_flow_mode"), &UnderGenLiquidFloodNode::get_flow_mode);
    ClassDB::bind_method(D_METHOD("set_flow_spread_limit", "limit"), &UnderGenLiquidFloodNode::set_flow_spread_limit);
    ClassDB::bind_method(D_METHOD("get_flow_spread_limit"), &UnderGenLiquidFloodNode::get_flow_spread_limit);
    ClassDB::bind_method(D_METHOD("set_source_spawn_chance", "chance"), &UnderGenLiquidFloodNode::set_source_spawn_chance);
    ClassDB::bind_method(D_METHOD("get_source_spawn_chance"), &UnderGenLiquidFloodNode::get_source_spawn_chance);
    ClassDB::bind_method(D_METHOD("set_max_sources", "max"), &UnderGenLiquidFloodNode::set_max_sources);
    ClassDB::bind_method(D_METHOD("get_max_sources"), &UnderGenLiquidFloodNode::get_max_sources);
    ClassDB::bind_method(D_METHOD("set_carve_basins", "carve"), &UnderGenLiquidFloodNode::set_carve_basins);
    ClassDB::bind_method(D_METHOD("get_carve_basins"), &UnderGenLiquidFloodNode::get_carve_basins);
    ClassDB::bind_method(D_METHOD("set_basin_radius", "radius"), &UnderGenLiquidFloodNode::set_basin_radius);
    ClassDB::bind_method(D_METHOD("get_basin_radius"), &UnderGenLiquidFloodNode::get_basin_radius);
    ClassDB::bind_method(D_METHOD("set_basin_depth", "depth"), &UnderGenLiquidFloodNode::set_basin_depth);
    ClassDB::bind_method(D_METHOD("get_basin_depth"), &UnderGenLiquidFloodNode::get_basin_depth);
    ClassDB::bind_method(D_METHOD("set_basin_carve_value", "value"), &UnderGenLiquidFloodNode::set_basin_carve_value);
    ClassDB::bind_method(D_METHOD("get_basin_carve_value"), &UnderGenLiquidFloodNode::get_basin_carve_value);

    ADD_PROPERTY(PropertyInfo(Variant::INT, "liquid_material_id", PROPERTY_HINT_RANGE, "0,255,1"), "set_liquid_material_id", "get_liquid_material_id");
    ADD_PROPERTY(PropertyInfo(Variant::DICTIONARY, "zone_flood_levels"), "set_zone_flood_levels", "get_zone_flood_levels");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "pool_radius"), "set_pool_radius", "get_pool_radius");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "pool_depth"), "set_pool_depth", "get_pool_depth");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "voxel_size"), "set_voxel_size", "get_voxel_size");

    ADD_PROPERTY(PropertyInfo(Variant::INT, "flow_mode", PROPERTY_HINT_ENUM, "Immediate Pool,Minecraft Sim"), "set_flow_mode", "get_flow_mode");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "flow_spread_limit", PROPERTY_HINT_RANGE, "1,15,1"), "set_flow_spread_limit", "get_flow_spread_limit");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "source_spawn_chance", PROPERTY_HINT_RANGE, "0.0,1.0,0.01"), "set_source_spawn_chance", "get_source_spawn_chance");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "max_sources"), "set_max_sources", "get_max_sources");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "carve_basins"), "set_carve_basins", "get_carve_basins");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "basin_radius"), "set_basin_radius", "get_basin_radius");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "basin_depth"), "set_basin_depth", "get_basin_depth");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "basin_carve_value"), "set_basin_carve_value", "get_basin_carve_value");

    BIND_ENUM_CONSTANT(FLOW_MODE_IMMEDIATE_POOL);
    BIND_ENUM_CONSTANT(FLOW_MODE_MINECRAFT_SIM);
}

void UnderGenLiquidFloodNode::set_liquid_material_id(int p_id) { liquid_material_id = Math::clamp(p_id, 0, 255); }
int UnderGenLiquidFloodNode::get_liquid_material_id() const { return liquid_material_id; }

void UnderGenLiquidFloodNode::set_zone_flood_levels(const Dictionary &p_levels) { zone_flood_levels = p_levels; }
Dictionary UnderGenLiquidFloodNode::get_zone_flood_levels() const { return zone_flood_levels; }

void UnderGenLiquidFloodNode::set_pool_radius(float p_radius) { pool_radius = p_radius > 0.0f ? p_radius : 0.0f; }
float UnderGenLiquidFloodNode::get_pool_radius() const { return pool_radius; }

void UnderGenLiquidFloodNode::set_pool_depth(float p_depth) { pool_depth = p_depth > 0.0f ? p_depth : 0.0f; }
float UnderGenLiquidFloodNode::get_pool_depth() const { return pool_depth; }

void UnderGenLiquidFloodNode::set_voxel_size(float p_size) { voxel_size = p_size > 0.0f ? p_size : 0.001f; }
float UnderGenLiquidFloodNode::get_voxel_size() const { return voxel_size; }

// Getters & Setters
void UnderGenLiquidFloodNode::set_flow_mode(FlowMode p_mode) { flow_mode = p_mode; }
UnderGenLiquidFloodNode::FlowMode UnderGenLiquidFloodNode::get_flow_mode() const { return flow_mode; }

void UnderGenLiquidFloodNode::set_flow_spread_limit(int p_limit) { flow_spread_limit = p_limit > 0 ? p_limit : 1; }
int UnderGenLiquidFloodNode::get_flow_spread_limit() const { return flow_spread_limit; }

void UnderGenLiquidFloodNode::set_source_spawn_chance(float p_chance) { source_spawn_chance = Math::clamp(p_chance, 0.0f, 1.0f); }
float UnderGenLiquidFloodNode::get_source_spawn_chance() const { return source_spawn_chance; }

void UnderGenLiquidFloodNode::set_max_sources(int p_max) { max_sources = p_max > 0 ? p_max : 1; }
int UnderGenLiquidFloodNode::get_max_sources() const { return max_sources; }

void UnderGenLiquidFloodNode::set_carve_basins(bool p_carve) { carve_basins = p_carve; }
bool UnderGenLiquidFloodNode::get_carve_basins() const { return carve_basins; }

void UnderGenLiquidFloodNode::set_basin_radius(float p_radius) { basin_radius = p_radius > 0.0f ? p_radius : 0.0f; }
float UnderGenLiquidFloodNode::get_basin_radius() const { return basin_radius; }

void UnderGenLiquidFloodNode::set_basin_depth(float p_depth) { basin_depth = p_depth > 0.0f ? p_depth : 0.0f; }
float UnderGenLiquidFloodNode::get_basin_depth() const { return basin_depth; }

void UnderGenLiquidFloodNode::set_basin_carve_value(float p_value) { basin_carve_value = p_value; }
float UnderGenLiquidFloodNode::get_basin_carve_value() const { return basin_carve_value; }

// Basin Carving Helper
void UnderGenLiquidFloodNode::_carve_basin(Ref<DensityGrid> grid, float v_size, const Vector3 &center, float radius, float depth, float surf_thresh, std::vector<Vector3i> &out_basin_cells) {
    if (grid.is_null()) return;

    Vector3i grid_center(
        (int)Math::round(center.x / v_size),
        (int)Math::round(center.y / v_size),
        (int)Math::round(center.z / v_size)
    );

    int r_cells = (int)Math::ceil(radius / v_size);
    int d_cells = (int)Math::ceil(depth / v_size);

    Vector3i dims = grid->get_grid_dimensions();
    int min_x = Math::max(0, grid_center.x - r_cells);
    int max_x = Math::min(dims.x - 1, grid_center.x + r_cells);
    int min_y = Math::max(0, grid_center.y - d_cells);
    int max_y = Math::min(dims.y - 1, grid_center.y);
    int min_z = Math::max(0, grid_center.z - r_cells);
    int max_z = Math::min(dims.z - 1, grid_center.z + r_cells);

    float r_sq = radius * radius;
    for (int z = min_z; z <= max_z; ++z) {
        float wz = z * v_size;
        float dz = wz - center.z;
        for (int x = min_x; x <= max_x; ++x) {
            float wx = x * v_size;
            float dx = wx - center.x;
            
            float h_dist_sq = dx * dx + dz * dz;
            if (h_dist_sq > r_sq) {
                continue;
            }
            
            float h_dist = Math::sqrt(h_dist_sq);
            float normalized_h = radius > 0.001f ? h_dist / radius : 0.0f;
            float bowl_depth_at_h = depth * (1.0f - normalized_h * normalized_h);
            float target_floor_y = center.y - bowl_depth_at_h;
            
            int end_y = Math::min(dims.y - 1, (int)Math::floor(center.y / v_size) - 1);
            int start_y = Math::max(0, (int)Math::floor(target_floor_y / v_size));
            if (start_y > end_y) {
                start_y = end_y;
            }
            
            for (int y = start_y; y <= end_y; ++y) {
                Vector3i pos(x, y, z);
                if (grid->get_cell(pos) > surf_thresh) {
                    grid->set_cell(pos, basin_carve_value);
                }
                out_basin_cells.push_back(pos);
            }
        }
    }
}

// Distance to nearest drop search
int UnderGenLiquidFloodNode::_find_distance_to_drop(Ref<DensityGrid> grid, float surf_thresh, const Vector3i &start, int max_dist) const {
    if (grid.is_null()) return 999;
    if (!grid->is_valid_position(start) || grid->get_cell(start) > surf_thresh) {
        return 999;
    }

    // Check if the cell directly below the start is air
    Vector3i down = start + Vector3i(0, -1, 0);
    if (grid->is_valid_position(down) && grid->get_cell(down) <= surf_thresh) {
        return 0;
    }

    if (max_dist <= 0) {
        return 999;
    }

    // Local horizontal BFS
    struct BFSNode {
        Vector3i pos;
        int dist;
    };
    std::queue<BFSNode> bfs_queue;
    std::vector<Vector3i> visited;

    bfs_queue.push({start, 0});
    visited.push_back(start);

    const Vector3i directions[4] = {
        Vector3i(1, 0, 0), Vector3i(-1, 0, 0),
        Vector3i(0, 0, 1), Vector3i(0, 0, -1)
    };

    while (!bfs_queue.empty()) {
        BFSNode curr = bfs_queue.front();
        bfs_queue.pop();

        if (curr.dist >= max_dist) {
            continue;
        }

        for (int i = 0; i < 4; ++i) {
            Vector3i next_pos = curr.pos + directions[i];

            // Check if visited
            bool is_visited = false;
            for (const auto &v : visited) {
                if (v == next_pos) {
                    is_visited = true;
                    break;
                }
            }
            if (is_visited) continue;

            if (grid->is_valid_position(next_pos) && grid->get_cell(next_pos) <= surf_thresh) {
                Vector3i next_down = next_pos + Vector3i(0, -1, 0);
                if (grid->is_valid_position(next_down) && grid->get_cell(next_down) <= surf_thresh) {
                    return curr.dist + 1;
                }

                visited.push_back(next_pos);
                bfs_queue.push({next_pos, curr.dist + 1});
            }
        }
    }

    return 999;
}

void UnderGenLiquidFloodNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    // Port 0: Input Generation Context
    Dictionary context = inputs.get(0, Dictionary());
    if (context.is_empty()) {
        UtilityFunctions::printerr("UnderGenLiquidFloodNode: Input context is empty.");
        return;
    }

    Ref<DensityGrid> grid = context.get("grid", Ref<DensityGrid>());
    if (grid.is_null()) {
        UtilityFunctions::printerr("UnderGenLiquidFloodNode: Grid not found in context.");
        context["flow_spread_limit"] = flow_spread_limit;
        outputs[0] = context;
        return;
    }

    // Retrieve voxel_size, looking up local property first, then context or default.
    float v_size = (voxel_size > 0.0f) ? voxel_size : (float)context.get("voxel_size", 1.0f);

    Vector3i dims = grid->get_grid_dimensions();
    int gsx = dims.x;
    int gsy = dims.y;
    int gsz = dims.z;
    float surf_thresh = grid->get_surface_threshold();

    // Map zone name strings to registered zone IDs.
    std::map<String, std::vector<int>> zone_name_to_ids;
    int zone_count = grid->get_zone_count();
    for (int zid = 0; zid < zone_count; ++zid) {
        String zone_name = grid->get_zone_name_by_id(zid);
        if (!zone_name.is_empty()) {
            zone_name_to_ids[zone_name].push_back(zid);
        }
    }

    // ────────────────────────────────────────────────────────────────────────
    // 1. Zone-Based Flooding
    // ────────────────────────────────────────────────────────────────────────
    if (!zone_flood_levels.is_empty()) {
        Array keys = zone_flood_levels.keys();
        for (int i = 0; i < keys.size(); ++i) {
            String zone_name = keys[i];
            Variant flood_config = zone_flood_levels[zone_name];

            float flood_height = 0.0f;
            bool relative_to_floor = false;

            if (flood_config.get_type() == Variant::DICTIONARY) {
                Dictionary config_dict = flood_config;
                flood_height = config_dict.get("height", 0.0f);
                relative_to_floor = config_dict.get("relative_to_floor", false);
            } else if (flood_config.get_type() == Variant::FLOAT || flood_config.get_type() == Variant::INT) {
                flood_height = (float)flood_config;
            } else {
                continue;
            }

            auto it = zone_name_to_ids.find(zone_name);
            if (it == zone_name_to_ids.end()) {
                continue; // Zone name not active in the current grid.
            }

            const std::vector<int> &zids = it->second;

            // Compute local floor height if relative_to_floor is requested.
            float floor_y = 99999.0f;
            if (relative_to_floor) {
                bool found_floor = false;
                for (int z = 0; z < gsz; ++z) {
                    for (int y = 0; y < gsy; ++y) {
                        for (int x = 0; x < gsx; ++x) {
                            Vector3i pos(x, y, z);
                            int cell_zid = grid->get_zone_at(pos);
                            bool match = false;
                            for (int zid : zids) {
                                  if (cell_zid == zid) {
                                      match = true;
                                      break;
                                  }
                            }
                            if (match && grid->get_cell(pos) <= surf_thresh) {
                                float world_y = y * v_size;
                                if (world_y < floor_y) {
                                    floor_y = world_y;
                                    found_floor = true;
                                }
                            }
                        }
                    }
                }
                if (!found_floor) {
                    floor_y = 0.0f;
                }
            }

            float final_flood_level_world = relative_to_floor ? (floor_y + flood_height) : flood_height;

            // Fill empty cells belonging to this zone below the flood level.
            for (int z = 0; z < gsz; ++z) {
                for (int y = 0; y < gsy; ++y) {
                    float cell_y_world = y * v_size;
                    if (cell_y_world > final_flood_level_world) {
                        continue;
                    }
                    for (int x = 0; x < gsx; ++x) {
                        Vector3i pos(x, y, z);
                        int cell_zid = grid->get_zone_at(pos);
                        bool match = false;
                        for (int zid : zids) {
                            if (cell_zid == zid) {
                                match = true;
                                break;
                            }
                        }
                        // Only fill air (carved density <= surface threshold).
                        if (match && grid->get_cell(pos) <= surf_thresh) {
                            grid->set_material_id(pos, 255); // Max level for flat zone flooding
                        }
                    }
                }
            }
        }
    }

    // ────────────────────────────────────────────────────────────────────────
    // 2. Point-Based Flooding (Local Pools or Minecraft Flow Sim)
    // ────────────────────────────────────────────────────────────────────────
    // Port 1: Optional Incoming PointSet
    Ref<UnderGenPointSet> in_points = inputs.get(1, Ref<UnderGenPointSet>());
    if (in_points.is_valid()) {
        if (flow_mode == FLOW_MODE_IMMEDIATE_POOL) {
            // Old Cylinder Pool Implementation
            if (pool_radius > 0.001f && pool_depth > 0.001f) {
                float r_sq = pool_radius * pool_radius;
                int point_count = in_points->get_point_count();

                for (int i = 0; i < point_count; ++i) {
                    Transform3D transform = in_points->get_point_transform(i);
                    Vector3 center = transform.origin;

                    int min_x = Math::max(0, (int)Math::floor((center.x - pool_radius) / v_size));
                    int max_x = Math::min(gsx - 1, (int)Math::ceil((center.x + pool_radius) / v_size));
                    int min_y = Math::max(0, (int)Math::floor((center.y - pool_depth) / v_size));
                    int max_y = Math::min(gsy - 1, (int)Math::ceil(center.y / v_size));
                    int min_z = Math::max(0, (int)Math::floor((center.z - pool_radius) / v_size));
                    int max_z = Math::min(gsz - 1, (int)Math::ceil((center.z + pool_radius) / v_size));

                    for (int z = min_z; z <= max_z; ++z) {
                        float wz = z * v_size;
                        float dz = wz - center.z;
                        for (int x = min_x; x <= max_x; ++x) {
                            float wx = x * v_size;
                            float dx = wx - center.x;
                            if (dx * dx + dz * dz > r_sq) {
                                continue;
                            }
                            for (int y = min_y; y <= max_y; ++y) {
                                Vector3i pos(x, y, z);
                                if (grid->get_cell(pos) <= surf_thresh) {
                                    grid->set_material_id(pos, 255); // Max level for legacy pool flooding
                                }
                            }
                        }
                    }
                }
            }
        } else if (flow_mode == FLOW_MODE_MINECRAFT_SIM) {
            // Seed RNG deterministically based on pipeline context
            int64_t seed = context.get("seed", (int64_t)12345);
            rng->set_seed(seed);

            // 2.1 Filter incoming points
            std::vector<Vector3> candidate_sources;
            int point_count = in_points->get_point_count();
            for (int i = 0; i < point_count; ++i) {
                Transform3D transform = in_points->get_point_transform(i);
                candidate_sources.push_back(transform.origin);
            }

            // Seeded shuffle of candidates to randomly distribute sources
            if (candidate_sources.size() > 1) {
                for (size_t i = candidate_sources.size() - 1; i > 0; --i) {
                    int j = rng->randi_range(0, (int)i);
                    std::swap(candidate_sources[i], candidate_sources[j]);
                }
            }

            // Select up to max_sources that pass spawn_chance check
            std::vector<Vector3> selected_sources;
            for (const Vector3 &pt : candidate_sources) {
                if (selected_sources.size() >= (size_t)max_sources) {
                    break;
                }
                if (rng->randf() <= source_spawn_chance) {
                    selected_sources.push_back(pt);
                }
            }

            // 2.2 Carve basins at source points
            std::vector<Vector3i> basin_cells;
            if (carve_basins && basin_radius > 0.001f && basin_depth > 0.001f) {
                for (const Vector3 &src_point : selected_sources) {
                    _carve_basin(grid, v_size, src_point, basin_radius, basin_depth, surf_thresh, basin_cells);
                }
            }

            // 2.3 Run Minecraft Flow Simulation
            struct WaterCell {
                Vector3i pos;
                int level;
            };
            std::queue<WaterCell> queue;
            std::vector<int> water_levels(gsx * gsy * gsz, 0);

            int max_level = flow_spread_limit + 1;

            // Push starting sources to the queue
            for (const Vector3 &src_point : selected_sources) {
                int sy = (int)Math::round(src_point.y / v_size);
                if (carve_basins && basin_radius > 0.001f && basin_depth > 0.001f) {
                    sy = (int)Math::floor(src_point.y / v_size) - 1;
                }
                Vector3i src_grid(
                    (int)Math::round(src_point.x / v_size),
                    sy,
                    (int)Math::round(src_point.z / v_size)
                );
                if (grid->is_valid_position(src_grid)) {
                    int idx = grid->get_index(src_grid);
                    if (water_levels[idx] < max_level) {
                        water_levels[idx] = max_level;
                        queue.push({src_grid, max_level});
                        grid->set_material_id(src_grid, 255); // Max level for starting sources
                    }
                }
            }

            // Fill all basin cells with water and add them to the queue at max level
            for (const Vector3i &cell : basin_cells) {
                if (grid->is_valid_position(cell)) {
                    int idx = grid->get_index(cell);
                    if (water_levels[idx] < max_level) {
                        water_levels[idx] = max_level;
                        queue.push({cell, max_level});
                        grid->set_material_id(cell, 255); // Max level inside basin
                    }
                }
            }

            int iterations = 0;
            const int max_iterations = 100000;

            const Vector3i directions[4] = {
                Vector3i(1, 0, 0), Vector3i(-1, 0, 0),
                Vector3i(0, 0, 1), Vector3i(0, 0, -1)
            };

            while (!queue.empty() && iterations < max_iterations) {
                iterations++;
                WaterCell curr = queue.front();
                queue.pop();

                // Only set material to liquid if the cell is air/carved
                if (grid->get_cell(curr.pos) <= surf_thresh) {
                    int mat_to_set = 255 - (max_level - curr.level);
                    grid->set_material_id(curr.pos, mat_to_set);
                }

                if (curr.level <= 0) {
                    continue;
                }

                // Check downward cell
                Vector3i down_pos = curr.pos + Vector3i(0, -1, 0);
                if (grid->is_valid_position(down_pos) && grid->get_cell(down_pos) <= surf_thresh) {
                    // Downward flow resets water level to maximum (waterfalls)
                    int down_idx = grid->get_index(down_pos);
                    if (water_levels[down_idx] < max_level) {
                        water_levels[down_idx] = max_level;
                        queue.push({down_pos, max_level});
                    }
                } else {
                    // Blocked below - spread horizontally
                    int dists[4];
                    int min_dist = 999;
                    for (int i = 0; i < 4; ++i) {
                        Vector3i next_pos = curr.pos + directions[i];
                        dists[i] = _find_distance_to_drop(grid, surf_thresh, next_pos, 5);
                        if (dists[i] < min_dist) {
                            min_dist = dists[i];
                        }
                    }

                    int next_level = curr.level - 1;
                    if (next_level > 0) {
                        for (int i = 0; i < 4; ++i) {
                            Vector3i next_pos = curr.pos + directions[i];

                            // Flow towards nearest drop if one was found
                            if (min_dist < 999 && dists[i] != min_dist) {
                                continue;
                            }

                            if (grid->is_valid_position(next_pos) && grid->get_cell(next_pos) <= surf_thresh) {
                                int next_idx = grid->get_index(next_pos);
                                if (water_levels[next_idx] < next_level) {
                                    water_levels[next_idx] = next_level;
                                    queue.push({next_pos, next_level});
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    context["flow_spread_limit"] = flow_spread_limit;
    outputs[0] = context;
}

} // namespace godot

