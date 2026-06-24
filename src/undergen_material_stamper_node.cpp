#include "undergen_material_stamper_node.h"
#include "density_grid.h"
#include <godot_cpp/variant/utility_functions.hpp>
#include <godot_cpp/core/math.hpp>
#include <map>

namespace godot {

UnderGenMaterialStamperNode::UnderGenMaterialStamperNode() {}
UnderGenMaterialStamperNode::~UnderGenMaterialStamperNode() {}

void UnderGenMaterialStamperNode::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_zone_material_map", "map"), &UnderGenMaterialStamperNode::set_zone_material_map);
    ClassDB::bind_method(D_METHOD("get_zone_material_map"), &UnderGenMaterialStamperNode::get_zone_material_map);
    ClassDB::bind_method(D_METHOD("set_default_material_id", "id"), &UnderGenMaterialStamperNode::set_default_material_id);
    ClassDB::bind_method(D_METHOD("get_default_material_id"), &UnderGenMaterialStamperNode::get_default_material_id);

    ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "zone_material_map", PROPERTY_HINT_TYPE_STRING,
                              vformat("%d/%d:%s", Variant::OBJECT, Variant::NIL, "ZoneMaterialEntry")),
                 "set_zone_material_map", "get_zone_material_map");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "default_material_id", PROPERTY_HINT_RANGE, "0,255,1"),
                 "set_default_material_id", "get_default_material_id");
}

void UnderGenMaterialStamperNode::set_zone_material_map(const TypedArray<ZoneMaterialEntry> &p_map) { zone_material_map = p_map; }
TypedArray<ZoneMaterialEntry> UnderGenMaterialStamperNode::get_zone_material_map() const { return zone_material_map; }
void UnderGenMaterialStamperNode::set_default_material_id(int p_id) { default_material_id = Math::clamp(p_id, 0, 255); }
int UnderGenMaterialStamperNode::get_default_material_id() const { return default_material_id; }

void UnderGenMaterialStamperNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    Dictionary context = inputs.get(0, Dictionary());
    if (context.is_empty()) {
        UtilityFunctions::printerr("UnderGenMaterialStamperNode: Input context is empty.");
        return;
    }

    Ref<DensityGrid> grid = context.get("grid", Ref<DensityGrid>());
    if (grid.is_null()) {
        UtilityFunctions::printerr("UnderGenMaterialStamperNode: Grid not found in context.");
        return;
    }

    int gsx = grid->get_grid_size_x();
    int gsy = grid->get_grid_size_y();
    int gsz = grid->get_grid_size_z();

    // ── Pass 1: Build zone name → zone id map ──────────────────────────
    std::map<String, int> zone_name_to_id;

    for (int z = 0; z < gsz; ++z) {
        for (int y = 0; y < gsy; ++y) {
            for (int x = 0; x < gsx; ++x) {
                int zid = grid->get_zone_at(Vector3i(x, y, z));
                if (zid > 0) {
                    String zname = grid->get_zone_name_by_id(zid);
                    if (zone_name_to_id.find(zname) == zone_name_to_id.end()) {
                        zone_name_to_id[zname] = zid;
                    }
                }
            }
        }
    }

    // ── Resolve zone_name → material_id entries, build zone_id → material lookup ──
    int max_zone_id = 0;
    for (const auto &pair : zone_name_to_id) {
        if (pair.second > max_zone_id) max_zone_id = pair.second;
    }

    std::vector<int> zone_id_to_material(max_zone_id + 1, default_material_id);

    for (int i = 0; i < zone_material_map.size(); ++i) {
        Ref<ZoneMaterialEntry> entry = zone_material_map[i];
        if (entry.is_null()) continue;

        String zone_name = entry->get_zone_name();
        int mat_id = entry->get_material_id();

        auto it = zone_name_to_id.find(zone_name);
        if (it != zone_name_to_id.end()) {
            int zid = it->second;
            if (zid >= 0 && zid <= max_zone_id) {
                zone_id_to_material[zid] = mat_id;
            }
        } else {
            UtilityFunctions::print("UnderGenMaterialStamperNode: Zone '", zone_name, "' not found in grid, skipping.");
        }
    }

    // ── Pass 2: Stamp material IDs onto every voxel ────────────────────
    int stamped = 0;
    for (int z = 0; z < gsz; ++z) {
        for (int y = 0; y < gsy; ++y) {
            for (int x = 0; x < gsx; ++x) {
                Vector3i pos(x, y, z);
                int zid = grid->get_zone_at(pos);
                if (zid >= 0 && zid < (int)zone_id_to_material.size()) {
                    grid->set_material_id(pos, zone_id_to_material[zid]);
                    ++stamped;
                }
            }
        }
    }

    UtilityFunctions::print("UnderGenMaterialStamperNode: Stamped ", stamped, " voxels.");
    outputs[0] = context;
}

} // namespace godot
