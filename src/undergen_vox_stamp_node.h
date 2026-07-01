#ifndef UNDERGEN_VOX_STAMP_NODE_H
#define UNDERGEN_VOX_STAMP_NODE_H

#include "undergen_node.h"
#include "level_gen_data.h"
#include "density_grid.h"
#include "ogt_vox.h"
#include <godot_cpp/variant/dictionary.hpp>
#include <map>
#include <vector>

namespace godot {

class UnderGenVoxStampNode : public UnderGenNode {
    GDCLASS(UnderGenVoxStampNode, UnderGenNode);

private:
    Dictionary vox_spawn_map;    // palette_index -> spawn_type String
    Dictionary vox_material_map; // palette_index -> material_id int
    bool vox_inverse_density = false;

    // Per-execution cache — cleared each time
    std::map<String, const ogt_vox_scene*> vox_cache;

    const ogt_vox_scene* _load_vox(const String &path);
    void _stamp_vox(DensityGrid* grid, const ResolvedRoom &room, const ogt_vox_scene* scene,
                    std::vector<Dictionary> &out_spawns);
    void _clear_vox_cache();

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

    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
};

} // namespace godot

#endif // UNDERGEN_VOX_STAMP_NODE_H
