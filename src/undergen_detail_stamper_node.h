#ifndef UNDERGEN_DETAIL_STAMPER_NODE_H
#define UNDERGEN_DETAIL_STAMPER_NODE_H

#include "undergen_node.h"
#include "undergen_point_set.h"
#include "ogt_vox.h"
#include "vox_spawn_entry.h"
#include "vox_material_entry.h"
#include <godot_cpp/classes/random_number_generator.hpp>
#include <godot_cpp/variant/dictionary.hpp>
#include <godot_cpp/variant/typed_array.hpp>
#include <map>
#include <vector>

namespace godot {

class UnderGenDetailStamperNode : public UnderGenNode {
    GDCLASS(UnderGenDetailStamperNode, UnderGenNode);

public:
    enum PivotPreset {
        PIVOT_BOTTOM_CENTER = 0,
        PIVOT_TOP_CENTER = 1,
        PIVOT_CENTER = 2,
        PIVOT_CUSTOM = 3
    };

private:
    String vox_path = "";
    float probability = 1.0f;
    bool align_with_normal = true;
    bool align_yaw_only = false;
    bool random_y_rotation = false;
    bool vox_inverse_density = false;
    PivotPreset pivot_preset = PIVOT_BOTTOM_CENTER;
    Vector3 custom_pivot = Vector3(0, 0, 0);
    int random_seed = 12345;

    Dictionary vox_spawn_map;
    Dictionary vox_material_map;

    TypedArray<VoxSpawnEntry> vox_spawn_entries;
    TypedArray<VoxMaterialEntry> vox_material_entries;

    Ref<RandomNumberGenerator> rng;

    // Caching loaded scenes
    std::map<String, const ogt_vox_scene*> vox_cache;

    const ogt_vox_scene* _load_vox(const String &path);
    void _clear_vox_cache();

    void _rebuild_maps_from_entries();
    void _on_entries_changed();

protected:
    static void _bind_methods();

public:
    UnderGenDetailStamperNode();
    virtual ~UnderGenDetailStamperNode();

    void set_vox_path(const String &p_path);
    String get_vox_path() const;

    void set_probability(float p_prob);
    float get_probability() const;

    void set_align_with_normal(bool p_align);
    bool get_align_with_normal() const;

    void set_align_yaw_only(bool p_enabled);
    bool get_align_yaw_only() const;

    void set_random_y_rotation(bool p_rand);
    bool get_random_y_rotation() const;

    void set_vox_inverse_density(bool p_enabled);
    bool get_vox_inverse_density() const;

    void set_pivot_preset(PivotPreset p_preset);
    PivotPreset get_pivot_preset() const;

    void set_custom_pivot(const Vector3 &p_pivot);
    Vector3 get_custom_pivot() const;

    void set_random_seed(int p_seed);
    int get_random_seed() const;

    void set_vox_spawn_entries(const TypedArray<VoxSpawnEntry> &p_entries);
    TypedArray<VoxSpawnEntry> get_vox_spawn_entries() const;
    void set_vox_material_entries(const TypedArray<VoxMaterialEntry> &p_entries);
    TypedArray<VoxMaterialEntry> get_vox_material_entries() const;

    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
};

} // namespace godot

VARIANT_ENUM_CAST(godot::UnderGenDetailStamperNode::PivotPreset);

#endif // UNDERGEN_DETAIL_STAMPER_NODE_H
