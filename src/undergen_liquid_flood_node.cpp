#include "undergen_liquid_flood_node.h"
#include "density_grid.h"
#include "undergen_point_set.h"
#include <godot_cpp/core/math.hpp>
#include <godot_cpp/variant/utility_functions.hpp>
#include <map>
#include <vector>

namespace godot {

UnderGenLiquidFloodNode::UnderGenLiquidFloodNode() {}
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

    ADD_PROPERTY(PropertyInfo(Variant::INT, "liquid_material_id", PROPERTY_HINT_RANGE, "0,255,1"), "set_liquid_material_id", "get_liquid_material_id");
    ADD_PROPERTY(PropertyInfo(Variant::DICTIONARY, "zone_flood_levels"), "set_zone_flood_levels", "get_zone_flood_levels");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "pool_radius"), "set_pool_radius", "get_pool_radius");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "pool_depth"), "set_pool_depth", "get_pool_depth");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "voxel_size"), "set_voxel_size", "get_voxel_size");
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
                            grid->set_material_id(pos, liquid_material_id);
                        }
                    }
                }
            }
        }
    }

    // ────────────────────────────────────────────────────────────────────────
    // 2. Point-Based Flooding (Local Pools)
    // ────────────────────────────────────────────────────────────────────────
    // Port 1: Optional Incoming PointSet
    Ref<UnderGenPointSet> in_points = inputs.get(1, Ref<UnderGenPointSet>());
    if (in_points.is_valid() && pool_radius > 0.001f && pool_depth > 0.001f) {
        float r_sq = pool_radius * pool_radius;
        int point_count = in_points->get_point_count();

        for (int i = 0; i < point_count; ++i) {
            Transform3D transform = in_points->get_point_transform(i);
            Vector3 center = transform.origin;

            // Convert world center and dimensions to grid bounds.
            int min_x = Math::max(0, (int)Math::floor((center.x - pool_radius) / v_size));
            int max_x = Math::min(gsx - 1, (int)Math::ceil((center.x + pool_radius) / v_size));
            int min_y = Math::max(0, (int)Math::floor((center.y - pool_depth) / v_size));
            int max_y = Math::min(gsy - 1, (int)Math::ceil(center.y / v_size)); // pool top is center Y
            int min_z = Math::max(0, (int)Math::floor((center.z - pool_radius) / v_size));
            int max_z = Math::min(gsz - 1, (int)Math::ceil((center.z + pool_radius) / v_size));

            for (int z = min_z; z <= max_z; ++z) {
                float wz = z * v_size;
                float dz = wz - center.z;
                for (int x = min_x; x <= max_x; ++x) {
                    float wx = x * v_size;
                    float dx = wx - center.x;
                    // Horizontal distance check.
                    if (dx * dx + dz * dz > r_sq) {
                        continue;
                    }
                    for (int y = min_y; y <= max_y; ++y) {
                        Vector3i pos(x, y, z);
                        if (grid->get_cell(pos) <= surf_thresh) {
                            grid->set_material_id(pos, liquid_material_id);
                        }
                    }
                }
            }
        }
    }

    outputs[0] = context;
}

} // namespace godot
