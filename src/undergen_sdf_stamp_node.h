#ifndef UNDERGEN_SDF_STAMP_NODE_H
#define UNDERGEN_SDF_STAMP_NODE_H

#include "undergen_node.h"
#include "level_gen_data.h"
#include "density_grid.h"
#include "vox_material_entry.h"
#include <godot_cpp/variant/dictionary.hpp>
#include <godot_cpp/variant/typed_array.hpp>
#include <godot_cpp/variant/quaternion.hpp>
#include <map>
#include <vector>

namespace godot {

struct SDFHermiteEdge {
    uint8_t edge_type; // 0: X, 1: Y, 2: Z
    uint16_t x, y, z;
    float t;
    Vector3 normal;
};

struct SDFMarker {
    String name;
    Vector3 pos;
    Quaternion rot;
};

struct SDFScene {
    uint32_t size_x = 0;
    uint32_t size_y = 0;
    uint32_t size_z = 0;
    float voxel_size = 1.0f;
    Vector3 pivot;
    std::vector<float> densities;
    std::vector<uint8_t> materials;
    std::vector<SDFHermiteEdge> hermite_edges;
    std::vector<SDFMarker> markers;
};

class UnderGenSdfStampNode : public UnderGenNode {
    GDCLASS(UnderGenSdfStampNode, UnderGenNode);

public:
    enum BlendMode {
        MODE_OVERWRITE = 0,
        MODE_MAX_UNION = 1,
        MODE_SUBTRACT_CARVE = 2
    };

private:
    Dictionary sdf_material_map; // sdf_material_id -> world_material_id
    int blend_mode = MODE_OVERWRITE;
    bool exclude_from_smoothing = false;
    bool exclude_from_warping = false;
    bool clear_room_air = true;

    TypedArray<VoxMaterialEntry> sdf_material_entries;

    // Per-execution cache — cleared each time
    std::map<String, SDFScene*> sdf_cache;

    SDFScene* _load_sdf(const String &path);
    void _clear_sdf_cache();

    void _rebuild_maps_from_entries();
    void _on_entries_changed();

protected:
    static void _bind_methods();

public:
    UnderGenSdfStampNode();
    virtual ~UnderGenSdfStampNode();

    static SDFScene* load_sdf_file(const String &path);
    static void stamp_sdf_scene(DensityGrid* grid, const ResolvedRoom &room, const SDFScene* scene,
                                const Dictionary &material_map, int blend_mode, int zone_id,
                                std::vector<Dictionary> &out_spawns, bool clear_room_air = true);

    void set_sdf_material_map(const Dictionary &p_map);
    Dictionary get_sdf_material_map() const;

    void set_blend_mode(int p_mode);
    int get_blend_mode() const;

    void set_exclude_from_smoothing(bool p_exclude);
    bool get_exclude_from_smoothing() const;

    void set_exclude_from_warping(bool p_exclude);
    bool get_exclude_from_warping() const;

    void set_clear_room_air(bool p_clear);
    bool get_clear_room_air() const;

    void set_sdf_material_entries(const TypedArray<VoxMaterialEntry> &p_entries);
    TypedArray<VoxMaterialEntry> get_sdf_material_entries() const;

    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
};

} // namespace godot

VARIANT_ENUM_CAST(UnderGenSdfStampNode::BlendMode);

#endif // UNDERGEN_SDF_STAMP_NODE_H
