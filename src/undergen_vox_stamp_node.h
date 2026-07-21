#ifndef UNDERGEN_VOX_STAMP_NODE_H
#define UNDERGEN_VOX_STAMP_NODE_H

#include "undergen_node.h"
#include "level_gen_data.h"
#include "density_grid.h"
#include "ogt_vox.h"
#include "vox_spawn_entry.h"
#include "vox_material_entry.h"
#include <godot_cpp/variant/dictionary.hpp>
#include <godot_cpp/variant/typed_array.hpp>
#include <map>
#include <vector>

namespace godot {

class UnderGenVoxStampNode : public UnderGenNode {
    GDCLASS(UnderGenVoxStampNode, UnderGenNode);

private:
    Dictionary vox_spawn_map;    // palette_index -> spawn_type String
    Dictionary vox_material_map; // palette_index -> material_id int
    bool vox_inverse_density = false;
    int connection_palette_index = -1;
    bool exclude_from_smoothing = false;
    bool exclude_from_warping = false;

    TypedArray<VoxSpawnEntry> vox_spawn_entries;
    TypedArray<VoxMaterialEntry> vox_material_entries;

    // Per-execution cache — cleared each time
    std::map<String, const ogt_vox_scene*> vox_cache;

    const ogt_vox_scene* _load_vox(const String &path);
    void _stamp_vox(DensityGrid* grid, const ResolvedRoom &room, const ogt_vox_scene* scene,
                    int zone_id, std::vector<Dictionary> &out_spawns, Array &out_connections);
    void _clear_vox_cache();

    void _rebuild_maps_from_entries();
    void _on_entries_changed();

protected:
    static void _bind_methods();

public:
    UnderGenVoxStampNode();
    virtual ~UnderGenVoxStampNode();

    void set_vox_spawn_map(const Dictionary &p_map);
    Dictionary get_vox_spawn_map() const;
    void set_vox_material_map(const Dictionary &p_map);
    Dictionary get_vox_material_map() const;
    void set_vox_inverse_density(bool p_enabled);
    bool get_vox_inverse_density() const;
    void set_connection_palette_index(int p_index);
    int get_connection_palette_index() const;
    void set_exclude_from_smoothing(bool p_exclude);
    bool get_exclude_from_smoothing() const;
    void set_exclude_from_warping(bool p_exclude);
    bool get_exclude_from_warping() const;

    void set_vox_spawn_entries(const TypedArray<VoxSpawnEntry> &p_entries);
    TypedArray<VoxSpawnEntry> get_vox_spawn_entries() const;
    void set_vox_material_entries(const TypedArray<VoxMaterialEntry> &p_entries);
    TypedArray<VoxMaterialEntry> get_vox_material_entries() const;

    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
};

} // namespace godot

#endif // UNDERGEN_VOX_STAMP_NODE_H
