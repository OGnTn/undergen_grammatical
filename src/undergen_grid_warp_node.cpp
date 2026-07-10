#include "undergen_grid_warp_node.h"
#include "density_grid.h"
#include "grid_parallel.h"
#include <godot_cpp/classes/resource.hpp>
#include <godot_cpp/variant/utility_functions.hpp>
#include <vector>
#include <cmath>
#include <cstring>

namespace godot {

inline float get_density_at(const std::vector<float> &data, int gsx, int gsy, int gsz, int x, int y, int z) {
    x = Math::clamp(x, 0, gsx - 1);
    y = Math::clamp(y, 0, gsy - 1);
    z = Math::clamp(z, 0, gsz - 1);
    return data[(z * gsy + y) * gsx + x];
}

inline float sample_trilinear(const std::vector<float> &data, int gsx, int gsy, int gsz, float x, float y, float z) {
    int x0 = (int)Math::floor(x);
    int y0 = (int)Math::floor(y);
    int z0 = (int)Math::floor(z);
    
    int x1 = x0 + 1;
    int y1 = y0 + 1;
    int z1 = z0 + 1;
    
    float tx = x - x0;
    float ty = y - y0;
    float tz = z - z0;
    
    float c000 = get_density_at(data, gsx, gsy, gsz, x0, y0, z0);
    float c100 = get_density_at(data, gsx, gsy, gsz, x1, y0, z0);
    float c010 = get_density_at(data, gsx, gsy, gsz, x0, y1, z0);
    float c110 = get_density_at(data, gsx, gsy, gsz, x1, y1, z0);
    
    float c001 = get_density_at(data, gsx, gsy, gsz, x0, y0, z1);
    float c101 = get_density_at(data, gsx, gsy, gsz, x1, y0, z1);
    float c011 = get_density_at(data, gsx, gsy, gsz, x0, y1, z1);
    float c111 = get_density_at(data, gsx, gsy, gsz, x1, y1, z1);
    
    float c00 = c000 * (1.0f - tx) + c100 * tx;
    float c10 = c010 * (1.0f - tx) + c110 * tx;
    float c01 = c001 * (1.0f - tx) + c101 * tx;
    float c11 = c011 * (1.0f - tx) + c111 * tx;
    
    float c0 = c00 * (1.0f - ty) + c10 * ty;
    float c1 = c01 * (1.0f - ty) + c11 * ty;
    
    return c0 * (1.0f - tz) + c1 * tz;
}

inline uint8_t sample_nearest_material(const std::vector<uint8_t> &data, int gsx, int gsy, int gsz, float x, float y, float z) {
    int rx = (int)Math::round(x);
    int ry = (int)Math::round(y);
    int rz = (int)Math::round(z);
    rx = Math::clamp(rx, 0, gsx - 1);
    ry = Math::clamp(ry, 0, gsy - 1);
    rz = Math::clamp(rz, 0, gsz - 1);
    return data[(rz * gsy + ry) * gsx + rx];
}

inline int sample_nearest_zone(const std::vector<int> &data, int gsx, int gsy, int gsz, float x, float y, float z) {
    int rx = (int)Math::round(x);
    int ry = (int)Math::round(y);
    int rz = (int)Math::round(z);
    rx = Math::clamp(rx, 0, gsx - 1);
    ry = Math::clamp(ry, 0, gsy - 1);
    rz = Math::clamp(rz, 0, gsz - 1);
    return data[(rz * gsy + ry) * gsx + rx];
}

UnderGenGridWarpNode::UnderGenGridWarpNode() {
    noise_generator.instantiate();
}

UnderGenGridWarpNode::~UnderGenGridWarpNode() {}

void UnderGenGridWarpNode::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_warp_amplitude", "amplitude"), &UnderGenGridWarpNode::set_warp_amplitude);
    ClassDB::bind_method(D_METHOD("get_warp_amplitude"), &UnderGenGridWarpNode::get_warp_amplitude);
    ClassDB::bind_method(D_METHOD("set_noise_frequency", "frequency"), &UnderGenGridWarpNode::set_noise_frequency);
    ClassDB::bind_method(D_METHOD("get_noise_frequency"), &UnderGenGridWarpNode::get_noise_frequency);
    ClassDB::bind_method(D_METHOD("set_noise_seed", "seed"), &UnderGenGridWarpNode::set_noise_seed);
    ClassDB::bind_method(D_METHOD("get_noise_seed"), &UnderGenGridWarpNode::get_noise_seed);
    ClassDB::bind_method(D_METHOD("set_noise_generator", "noise"), &UnderGenGridWarpNode::set_noise_generator);
    ClassDB::bind_method(D_METHOD("get_noise_generator"), &UnderGenGridWarpNode::get_noise_generator);

    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "warp_amplitude"), "set_warp_amplitude", "get_warp_amplitude");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "noise_frequency"), "set_noise_frequency", "get_noise_frequency");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "noise_seed"), "set_noise_seed", "get_noise_seed");
    ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "noise_generator", PROPERTY_HINT_RESOURCE_TYPE, "FastNoiseLite"), "set_noise_generator", "get_noise_generator");
}

void UnderGenGridWarpNode::set_warp_amplitude(float p_amplitude) { warp_amplitude = p_amplitude; }
float UnderGenGridWarpNode::get_warp_amplitude() const { return warp_amplitude; }
void UnderGenGridWarpNode::set_noise_frequency(float p_freq) { noise_frequency = p_freq; }
float UnderGenGridWarpNode::get_noise_frequency() const { return noise_frequency; }
void UnderGenGridWarpNode::set_noise_seed(int p_seed) { noise_seed = p_seed; }
int UnderGenGridWarpNode::get_noise_seed() const { return noise_seed; }
void UnderGenGridWarpNode::set_noise_generator(const Ref<FastNoiseLite> &p_noise) { noise_generator = p_noise; }
Ref<FastNoiseLite> UnderGenGridWarpNode::get_noise_generator() const { return noise_generator; }

void UnderGenGridWarpNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    Dictionary context = inputs.get(0, Dictionary());
    if (context.is_empty()) return;

    Ref<DensityGrid> grid = context.get("grid", Ref<DensityGrid>());
    if (grid.is_null() || !noise_generator.is_valid()) {
        outputs[0] = context;
        return;
    }

    noise_generator->set_seed(noise_seed);
    noise_generator->set_frequency(noise_frequency);

    int gsx = grid->get_grid_size_x();
    int gsy = grid->get_grid_size_y();
    int gsz = grid->get_grid_size_z();
    int64_t total_size = grid->get_total_cell_count();
    if (total_size <= 0) {
        outputs[0] = context;
        return;
    }

    PackedFloat32Array &density_data = grid->get_density_data_rw();
    PackedByteArray &material_data = grid->get_material_data_rw();
    PackedInt32Array &zone_data = grid->get_zone_data_rw();

    if (density_data.size() < total_size || material_data.size() < total_size || zone_data.size() < total_size) {
        outputs[0] = context;
        return;
    }

    Array rooms_arr = context.get("rooms", Array());
    struct RoomBounds {
        Vector3i min;
        Vector3i max;
    };
    std::vector<RoomBounds> excluded_rooms;

    for (int i = 0; i < rooms_arr.size(); ++i) {
        Dictionary r_dict = rooms_arr[i];
        if (r_dict.get("exclude_from_warping", false)) {
            Vector3i pos = r_dict.get("position", Vector3i());
            Vector3i size = r_dict.get("size", Vector3i());
            RoomBounds b;
            b.min = pos;
            b.max = pos + size;
            excluded_rooms.push_back(b);
        }
    }

    // Clone original grids to read from
    std::vector<float> orig_density(total_size);
    std::memcpy(orig_density.data(), density_data.ptr(), total_size * sizeof(float));

    std::vector<uint8_t> orig_material(total_size);
    std::memcpy(orig_material.data(), material_data.ptr(), total_size * sizeof(uint8_t));

    std::vector<int> orig_zone(total_size);
    std::memcpy(orig_zone.data(), zone_data.ptr(), total_size * sizeof(int));

    float *density_rw = density_data.ptrw();
    uint8_t *material_rw = material_data.ptrw();
    int *zone_rw = zone_data.ptrw();

    int worker_count = grid_parallel_worker_count(gsz, total_size);
    std::vector<Ref<FastNoiseLite>> worker_noise(worker_count);
    for (int i = 0; i < worker_count; ++i) {
        Ref<Resource> duplicated = noise_generator->duplicate(true);
        worker_noise[i] = duplicated;
        if (worker_noise[i].is_null()) {
            worker_noise[i].instantiate();
        }
        worker_noise[i]->set_seed(noise_seed);
        worker_noise[i]->set_frequency(noise_frequency);
    }

    parallel_for_z(gsz, total_size, [&](int worker_index, int z_begin, int z_end) {
        Ref<FastNoiseLite> local_noise = worker_noise[worker_index];
        for (int z = z_begin; z < z_end; ++z) {
            int slice_offset = z * gsy * gsx;
            for (int y = 0; y < gsy; ++y) {
                int row_offset = slice_offset + y * gsx;
                for (int x = 0; x < gsx; ++x) {
                    int idx = row_offset + x;

                    // Check if cell is in any excluded room
                    bool is_excluded = false;
                    for (const auto& r : excluded_rooms) {
                        if (x >= r.min.x && x < r.max.x &&
                            y >= r.min.y && y < r.max.y &&
                            z >= r.min.z && z < r.max.z) {
                            is_excluded = true;
                            break;
                        }
                    }

                    if (is_excluded) {
                        density_rw[idx] = orig_density[idx];
                        material_rw[idx] = orig_material[idx];
                        zone_rw[idx] = orig_zone[idx];
                    } else {
                        // Sample 3 independent channels of noise to warp the 3 coordinate axes
                        float dx = local_noise->get_noise_3d((float)x, (float)y, (float)z) * warp_amplitude;
                        float dy = local_noise->get_noise_3d((float)x + 1000.0f, (float)y + 1000.0f, (float)z + 1000.0f) * warp_amplitude;
                        float dz = local_noise->get_noise_3d((float)x - 1000.0f, (float)y - 1000.0f, (float)z - 1000.0f) * warp_amplitude;

                        float warped_x = (float)x + dx;
                        float warped_y = (float)y + dy;
                        float warped_z = (float)z + dz;

                        // Write warped values back
                        density_rw[idx] = sample_trilinear(orig_density, gsx, gsy, gsz, warped_x, warped_y, warped_z);
                        material_rw[idx] = sample_nearest_material(orig_material, gsx, gsy, gsz, warped_x, warped_y, warped_z);
                        zone_rw[idx] = sample_nearest_zone(orig_zone, gsx, gsy, gsz, warped_x, warped_y, warped_z);
                    }
                }
            }
        }
    });

    outputs[0] = context;
}

} // namespace godot
