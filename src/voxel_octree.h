#ifndef VOXEL_OCTREE_H
#define VOXEL_OCTREE_H

#include "svo_node.h"
#include <godot_cpp/variant/vector3i.hpp>
#include <godot_cpp/variant/vector3.hpp>

namespace godot {

class VoxelOctree {
private:
    SVONode* root = nullptr;
    int root_size = 16;
    int max_depth = 4; // 16 -> 8 -> 4 -> 2 -> 1 voxel size

    SVONode* _get_or_create_leaf(const Vector3i &pos, int &out_local_x, int &out_local_y, int &out_local_z);
    const SVONode* _find_leaf(const Vector3i &pos, int &out_local_x, int &out_local_y, int &out_local_z) const;
    void _try_collapse(SVONode* node);

public:
    VoxelOctree(int p_size = 16, float default_density = 1.0f, uint8_t default_material = 0, int32_t default_zone = 0);
    ~VoxelOctree();

    void initialize(int p_size, float default_density = 1.0f, uint8_t default_material = 0, int32_t default_zone = 0);

    float get_density(const Vector3i &pos) const;
    void set_density(const Vector3i &pos, float val);

    uint8_t get_material(const Vector3i &pos) const;
    void set_material(const Vector3i &pos, uint8_t mat_id);

    int32_t get_zone(const Vector3i &pos) const;
    void set_zone(const Vector3i &pos, int32_t zone_id);

    bool is_uniform() const { return root ? root->is_homogeneous() : true; }
    SVONode* get_root() { return root; }
    const SVONode* get_root() const { return root; }

    int get_root_size() const { return root_size; }
    int get_max_depth() const { return max_depth; }

    // Fast hierarchical raycast
    bool raycast(const Vector3 &origin, const Vector3 &dir, float max_dist, Vector3 &out_hit_pos, Vector3 &out_hit_normal, float iso_surface = 0.0f) const;
};

} // namespace godot

#endif // VOXEL_OCTREE_H
