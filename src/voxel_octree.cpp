#include "voxel_octree.h"
#include <algorithm>
#include <cmath>

namespace godot {

VoxelOctree::VoxelOctree(int p_size, float default_density, uint8_t default_material, int32_t default_zone) {
    initialize(p_size, default_density, default_material, default_zone);
}

VoxelOctree::~VoxelOctree() {
    if (root) {
        delete root;
        root = nullptr;
    }
}

void VoxelOctree::initialize(int p_size, float default_density, uint8_t default_material, int32_t default_zone) {
    if (root) {
        delete root;
    }
    root_size = p_size;
    max_depth = 0;
    int cur = p_size;
    while (cur > 1) {
        cur >>= 1;
        max_depth++;
    }
    root = new SVONode(default_density, default_material, default_zone);
}

const SVONode* VoxelOctree::_find_leaf(const Vector3i &pos, int &out_local_x, int &out_local_y, int &out_local_z) const {
    if (!root) return nullptr;
    if (pos.x < 0 || pos.x >= root_size || pos.y < 0 || pos.y >= root_size || pos.z < 0 || pos.z >= root_size) {
        return nullptr;
    }

    const SVONode* curr = root;
    int curr_size = root_size;
    Vector3i curr_min(0, 0, 0);

    while (curr && curr->type == SVO_NODE_BRANCH) {
        curr_size >>= 1;
        int child_x = (pos.x >= curr_min.x + curr_size) ? 1 : 0;
        int child_y = (pos.y >= curr_min.y + curr_size) ? 1 : 0;
        int child_z = (pos.z >= curr_min.z + curr_size) ? 1 : 0;

        if (child_x) curr_min.x += curr_size;
        if (child_y) curr_min.y += curr_size;
        if (child_z) curr_min.z += curr_size;

        int child_idx = child_x | (child_y << 1) | (child_z << 2);
        curr = curr->children[child_idx];
    }

    out_local_x = pos.x - curr_min.x;
    out_local_y = pos.y - curr_min.y;
    out_local_z = pos.z - curr_min.z;
    return curr;
}

SVONode* VoxelOctree::_get_or_create_leaf(const Vector3i &pos, int &out_local_x, int &out_local_y, int &out_local_z) {
    if (!root) {
        initialize(root_size);
    }
    if (pos.x < 0 || pos.x >= root_size || pos.y < 0 || pos.y >= root_size || pos.z < 0 || pos.z >= root_size) {
        return nullptr;
    }

    SVONode* curr = root;
    int curr_size = root_size;
    Vector3i curr_min(0, 0, 0);

    while (curr_size > 1) {
        if (curr->type == SVO_NODE_UNIFORM) {
            // Split uniform node into branch with 8 uniform children
            float d = curr->uniform_density;
            uint8_t m = curr->uniform_material;
            int32_t z = curr->uniform_zone;

            curr->type = SVO_NODE_BRANCH;
            for (int i = 0; i < 8; ++i) {
                curr->children[i] = new SVONode(d, m, z);
            }
        } else if (curr->type == SVO_NODE_LEAF_GRID) {
            // Reached leaf grid block
            break;
        }

        curr_size >>= 1;
        int child_x = (pos.x >= curr_min.x + curr_size) ? 1 : 0;
        int child_y = (pos.y >= curr_min.y + curr_size) ? 1 : 0;
        int child_z = (pos.z >= curr_min.z + curr_size) ? 1 : 0;

        if (child_x) curr_min.x += curr_size;
        if (child_y) curr_min.y += curr_size;
        if (child_z) curr_min.z += curr_size;

        int child_idx = child_x | (child_y << 1) | (child_z << 2);
        if (!curr->children[child_idx]) {
            curr->children[child_idx] = new SVONode(curr->uniform_density, curr->uniform_material, curr->uniform_zone);
        }
        curr = curr->children[child_idx];
    }

    out_local_x = pos.x - curr_min.x;
    out_local_y = pos.y - curr_min.y;
    out_local_z = pos.z - curr_min.z;
    return curr;
}

float VoxelOctree::get_density(const Vector3i &pos) const {
    int lx = 0, ly = 0, lz = 0;
    const SVONode* leaf = _find_leaf(pos, lx, ly, lz);
    if (!leaf) return 1.0f;

    if (leaf->type == SVO_NODE_LEAF_GRID && !leaf->density_block.empty()) {
        int idx = lx + 16 * (ly + 16 * lz); // Assuming max block or local index
        if (idx >= 0 && idx < (int)leaf->density_block.size()) {
            return leaf->density_block[idx];
        }
    }
    return leaf->uniform_density;
}

void VoxelOctree::set_density(const Vector3i &pos, float val) {
    int lx = 0, ly = 0, lz = 0;
    SVONode* leaf = _get_or_create_leaf(pos, lx, ly, lz);
    if (!leaf) return;

    if (leaf->type == SVO_NODE_LEAF_GRID) {
        int idx = lx + 16 * (ly + 16 * lz);
        if (idx >= 0 && idx < (int)leaf->density_block.size()) {
            leaf->density_block[idx] = val;
        }
    } else {
        leaf->uniform_density = val;
    }
}

uint8_t VoxelOctree::get_material(const Vector3i &pos) const {
    int lx = 0, ly = 0, lz = 0;
    const SVONode* leaf = _find_leaf(pos, lx, ly, lz);
    if (!leaf) return 0;

    if (leaf->type == SVO_NODE_LEAF_GRID && !leaf->material_block.empty()) {
        int idx = lx + 16 * (ly + 16 * lz);
        if (idx >= 0 && idx < (int)leaf->material_block.size()) {
            return leaf->material_block[idx];
        }
    }
    return leaf->uniform_material;
}

void VoxelOctree::set_material(const Vector3i &pos, uint8_t mat_id) {
    int lx = 0, ly = 0, lz = 0;
    SVONode* leaf = _get_or_create_leaf(pos, lx, ly, lz);
    if (!leaf) return;

    if (leaf->type == SVO_NODE_LEAF_GRID) {
        int idx = lx + 16 * (ly + 16 * lz);
        if (idx >= 0 && idx < (int)leaf->material_block.size()) {
            leaf->material_block[idx] = mat_id;
        }
    } else {
        leaf->uniform_material = mat_id;
    }
}

int32_t VoxelOctree::get_zone(const Vector3i &pos) const {
    int lx = 0, ly = 0, lz = 0;
    const SVONode* leaf = _find_leaf(pos, lx, ly, lz);
    if (!leaf) return 0;

    if (leaf->type == SVO_NODE_LEAF_GRID && !leaf->zone_block.empty()) {
        int idx = lx + 16 * (ly + 16 * lz);
        if (idx >= 0 && idx < (int)leaf->zone_block.size()) {
            return leaf->zone_block[idx];
        }
    }
    return leaf->uniform_zone;
}

void VoxelOctree::set_zone(const Vector3i &pos, int32_t zone_id) {
    int lx = 0, ly = 0, lz = 0;
    SVONode* leaf = _get_or_create_leaf(pos, lx, ly, lz);
    if (!leaf) return;

    if (leaf->type == SVO_NODE_LEAF_GRID) {
        int idx = lx + 16 * (ly + 16 * lz);
        if (idx >= 0 && idx < (int)leaf->zone_block.size()) {
            leaf->zone_block[idx] = zone_id;
        }
    } else {
        leaf->uniform_zone = zone_id;
    }
}

bool VoxelOctree::raycast(const Vector3 &origin, const Vector3 &dir, float max_dist, Vector3 &out_hit_pos, Vector3 &out_hit_normal, float iso_surface) const {
    Vector3 norm_dir = dir.normalized();
    float step_size = 0.5f;
    float traveled = 0.0f;
    Vector3 curr_pos = origin;

    while (traveled <= max_dist) {
        Vector3i grid_p((int)std::floor(curr_pos.x), (int)std::floor(curr_pos.y), (int)std::floor(curr_pos.z));
        float den = get_density(grid_p);

        if (den <= iso_surface) {
            out_hit_pos = curr_pos;
            // Approximate surface normal using central difference
            float dx = get_density(grid_p + Vector3i(1, 0, 0)) - get_density(grid_p - Vector3i(1, 0, 0));
            float dy = get_density(grid_p + Vector3i(0, 1, 0)) - get_density(grid_p - Vector3i(0, 1, 0));
            float dz = get_density(grid_p + Vector3i(0, 0, 1)) - get_density(grid_p - Vector3i(0, 0, 1));
            Vector3 norm(dx, dy, dz);
            out_hit_normal = norm.length_squared() > 1e-6f ? norm.normalized() : -norm_dir;
            return true;
        }

        curr_pos += norm_dir * step_size;
        traveled += step_size;
    }
    return false;
}

} // namespace godot
