#include "undergen_material_stamper_node.h"
#include "density_grid.h"
#include "grid_parallel.h"

#include <godot_cpp/core/math.hpp>
#include <godot_cpp/variant/utility_functions.hpp>

#include <atomic>
#include <map>
#include <vector>

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
    int64_t total_size = grid->get_total_cell_count();
    if (total_size <= 0 ||
        grid->get_zone_data().size() < total_size ||
        grid->get_material_data().size() < total_size) {
        outputs[0] = context;
        return;
    }

    std::map<String, std::vector<int>> zone_name_to_ids;
    int zone_count = grid->get_zone_count();
    for (int zid = 1; zid < zone_count; ++zid) {
        String zone_name = grid->get_zone_name_by_id(zid);
        if (!zone_name.is_empty()) {
            zone_name_to_ids[zone_name].push_back(zid);
        }
    }

    int total_ids = 0;
    for (const auto &pair : zone_name_to_ids) {
        total_ids += (int)pair.second.size();
    }
    UtilityFunctions::print("UnderGenMaterialStamperNode: Zones registered (", (int)zone_name_to_ids.size(), " names, ", total_ids, " ids):");
    for (const auto &pair : zone_name_to_ids) {
        String ids_str;
        for (size_t j = 0; j < pair.second.size(); ++j) {
            if (j > 0) ids_str += ",";
            ids_str += String::num_int64(pair.second[j]);
        }
        UtilityFunctions::print("   name=\"", pair.first, "\"  ids=[", ids_str, "]");
    }

    int max_zone_id = Math::max(0, zone_count - 1);
    std::vector<int> zone_id_to_material(max_zone_id + 1, default_material_id);

    for (int i = 0; i < zone_material_map.size(); ++i) {
        Ref<ZoneMaterialEntry> entry = zone_material_map[i];
        if (entry.is_null()) continue;

        String zone_name = entry->get_zone_name();
        int mat_id = entry->get_material_id();

        auto it = zone_name_to_ids.find(zone_name);
        if (it != zone_name_to_ids.end()) {
            for (int zid : it->second) {
                if (zid >= 0 && zid <= max_zone_id) {
                    zone_id_to_material[zid] = mat_id;
                }
            }
            UtilityFunctions::print("UnderGenMaterialStamperNode: Mapped zone \"", zone_name, "\" (", (int)it->second.size(), " ids) -> material ", mat_id);
        } else {
            UtilityFunctions::print("UnderGenMaterialStamperNode: Zone '", zone_name, "' not registered, skipping.");
        }
    }

    const PackedInt32Array &zone_array = grid->get_zone_data();
    PackedByteArray &material_array = grid->get_material_data_rw();
    const int32_t *zone_data = zone_array.ptr();
    uint8_t *material_data = material_array.ptrw();
    std::atomic<int64_t> stamped(0);

    parallel_for_z(gsz, total_size, [&](int, int z_begin, int z_end) {
        int64_t local_stamped = 0;
        for (int z = z_begin; z < z_end; ++z) {
            int slice_offset = z * gsy * gsx;
            for (int y = 0; y < gsy; ++y) {
                int row_offset = slice_offset + y * gsx;
                for (int x = 0; x < gsx; ++x) {
                    int idx = row_offset + x;
                    int zid = zone_data[idx];
                    if (zid >= 0 && zid < (int)zone_id_to_material.size()) {
                        material_data[idx] = (uint8_t)Math::clamp(zone_id_to_material[zid], 0, 255);
                        ++local_stamped;
                    }
                }
            }
        }
        stamped.fetch_add(local_stamped, std::memory_order_relaxed);
    });

    UtilityFunctions::print("UnderGenMaterialStamperNode: Stamped ", (int64_t)stamped.load(), " voxels.");
    outputs[0] = context;
}

} // namespace godot
