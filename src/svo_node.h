#ifndef SVO_NODE_H
#define SVO_NODE_H

#include <cstdint>
#include <vector>
#include <godot_cpp/variant/vector3i.hpp>
#include <godot_cpp/variant/vector3.hpp>

namespace godot {

enum SVONodeType : uint8_t {
    SVO_NODE_UNIFORM = 0,
    SVO_NODE_BRANCH = 1,
    SVO_NODE_LEAF_GRID = 2
};

struct SVONode {
    SVONodeType type = SVO_NODE_UNIFORM;
    float uniform_density = 1.0f;
    uint8_t uniform_material = 0;
    int32_t uniform_zone = 0;

    // Children for SVO_NODE_BRANCH (8 octants)
    SVONode* children[8] = { nullptr };

    // Dense storage if node is SVO_NODE_LEAF_GRID (e.g., 8x8x8 or 16x16x16 block)
    std::vector<float> density_block;
    std::vector<uint8_t> material_block;
    std::vector<int32_t> zone_block;

    SVONode() {
        for (int i = 0; i < 8; ++i) {
            children[i] = nullptr;
        }
    }

    SVONode(float p_density, uint8_t p_material = 0, int32_t p_zone = 0)
        : type(SVO_NODE_UNIFORM), uniform_density(p_density), uniform_material(p_material), uniform_zone(p_zone) {
        for (int i = 0; i < 8; ++i) {
            children[i] = nullptr;
        }
    }

    ~SVONode() {
        clear_children();
    }

    void clear_children() {
        for (int i = 0; i < 8; ++i) {
            if (children[i]) {
                delete children[i];
                children[i] = nullptr;
            }
        }
    }

    bool is_homogeneous() const {
        return type == SVO_NODE_UNIFORM;
    }
};

} // namespace godot

#endif // SVO_NODE_H
