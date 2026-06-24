#ifndef UNDERGEN_MATERIAL_STAMPER_NODE_H
#define UNDERGEN_MATERIAL_STAMPER_NODE_H

#include "undergen_node.h"
#include "zone_material_entry.h"
#include <godot_cpp/variant/typed_array.hpp>

namespace godot {

// Stamps material IDs into the density grid based on zone-to-material mapping.
// Each ZoneMaterialEntry pairs a zone name with a material ID;
// zones not explicitly listed receive the default_material_id.
class UnderGenMaterialStamperNode : public UnderGenNode {
    GDCLASS(UnderGenMaterialStamperNode, UnderGenNode);

private:
    TypedArray<ZoneMaterialEntry> zone_material_map;

    // Material ID to use for voxels whose zone is not in the map
    int default_material_id = 0;

protected:
    static void _bind_methods();

public:
    UnderGenMaterialStamperNode();
    virtual ~UnderGenMaterialStamperNode();

    void set_zone_material_map(const TypedArray<ZoneMaterialEntry> &p_map);
    TypedArray<ZoneMaterialEntry> get_zone_material_map() const;
    void set_default_material_id(int p_id);
    int get_default_material_id() const;

    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
};

} // namespace godot

#endif // UNDERGEN_MATERIAL_STAMPER_NODE_H
