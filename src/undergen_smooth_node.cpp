#include "undergen_smooth_node.h"
#include "density_grid.h"
#include <godot_cpp/variant/utility_functions.hpp>
#include <vector>
#include <cstring>

namespace godot {

UnderGenSmoothNode::UnderGenSmoothNode() {}
UnderGenSmoothNode::~UnderGenSmoothNode() {}

void UnderGenSmoothNode::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_smoothing_strength", "strength"), &UnderGenSmoothNode::set_smoothing_strength);
    ClassDB::bind_method(D_METHOD("get_smoothing_strength"), &UnderGenSmoothNode::get_smoothing_strength);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "smoothing_strength"), "set_smoothing_strength", "get_smoothing_strength");
}

void UnderGenSmoothNode::set_smoothing_strength(int p_strength) { smoothing_strength = p_strength; }
int UnderGenSmoothNode::get_smoothing_strength() const { return smoothing_strength; }

void UnderGenSmoothNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    Dictionary context = inputs.get(0, Dictionary());
    if (context.is_empty()) return;

    Ref<DensityGrid> grid = context.get("grid", Ref<DensityGrid>());
    if (grid.is_null() || smoothing_strength <= 0) {
        outputs[0] = context;
        return;
    }

    int gsx = grid->get_grid_size_x();
    int gsy = grid->get_grid_size_y();
    int gsz = grid->get_grid_size_z();
    size_t total_size = (size_t)gsx * gsy * gsz;
    if (total_size == 0) { outputs[0] = context; return; }

    PackedFloat32Array data = grid->get_world_density_grid();
    float* grid_data = data.ptrw();
    std::vector<float> buffer(total_size);
    float* temp_data = buffer.data();

    int R = smoothing_strength;

    // Pass 1: Blur X
    for (int z = 0; z < gsz; ++z) {
        for (int y = 0; y < gsy; ++y) {
            int row_offset = (z * gsy + y) * gsx;
            for (int x = 0; x < gsx; ++x) {
                float sum = 0.0f;
                int count = 0;
                for (int k = -R; k <= R; ++k) {
                    int nx = x + k;
                    sum += (nx >= 0 && nx < gsx) ? grid_data[row_offset + nx] : WORLD_SOLID_VALUE;
                    count++;
                }
                temp_data[row_offset + x] = (count > 0) ? sum / count : grid_data[row_offset + x];
            }
        }
    }

    // Pass 2: Blur Y
    for (int z = 0; z < gsz; ++z) {
        int slice_offset = z * gsy * gsx;
        for (int x = 0; x < gsx; ++x) {
            for (int y = 0; y < gsy; ++y) {
                float sum = 0.0f;
                int count = 0;
                for (int k = -R; k <= R; ++k) {
                    int ny = y + k;
                    sum += (ny >= 0 && ny < gsy) ? temp_data[slice_offset + ny * gsx + x] : WORLD_SOLID_VALUE;
                    count++;
                }
                grid_data[slice_offset + y * gsx + x] = (count > 0) ? sum / count : temp_data[slice_offset + y * gsx + x];
            }
        }
    }

    // Pass 3: Blur Z
    int stride_z = gsx * gsy;
    for (int y = 0; y < gsy; ++y) {
        for (int x = 0; x < gsx; ++x) {
            int col_offset = y * gsx + x;
            for (int z = 0; z < gsz; ++z) {
                float sum = 0.0f;
                int count = 0;
                for (int k = -R; k <= R; ++k) {
                    int nz = z + k;
                    sum += (nz >= 0 && nz < gsz) ? grid_data[nz * stride_z + col_offset] : WORLD_SOLID_VALUE;
                    count++;
                }
                temp_data[z * stride_z + col_offset] = (count > 0) ? sum / count : grid_data[z * stride_z + col_offset];
            }
        }
    }
    // Enforce solid boundary casing to prevent holes in the terrain after smoothing
    for (int z = 0; z < gsz; ++z) {
        for (int y = 0; y < gsy; ++y) {
            temp_data[(z * gsy + y) * gsx + 0] = WORLD_SOLID_VALUE;
            temp_data[(z * gsy + y) * gsx + (gsx - 1)] = WORLD_SOLID_VALUE;
        }
    }
    for (int z = 0; z < gsz; ++z) {
        for (int x = 0; x < gsx; ++x) {
            temp_data[(z * gsy + 0) * gsx + x] = WORLD_SOLID_VALUE;
            temp_data[(z * gsy + (gsy - 1)) * gsx + x] = WORLD_SOLID_VALUE;
        }
    }
    for (int y = 0; y < gsy; ++y) {
        for (int x = 0; x < gsx; ++x) {
            temp_data[(0 * gsy + y) * gsx + x] = WORLD_SOLID_VALUE;
            temp_data[((gsz - 1) * gsy + y) * gsx + x] = WORLD_SOLID_VALUE;
        }
    }

    memcpy(grid_data, temp_data, total_size * sizeof(float));

    // Write back
    grid->set_world_density_grid(data);

    outputs[0] = context;
}

} // namespace godot
