#include "liquid_generator.h"
#include <godot_cpp/variant/utility_functions.hpp>
#include <godot_cpp/core/math.hpp>
#include <queue>
#include <unordered_set>
#include <vector>

namespace godot {

Ref<DensityGrid> LiquidGenerator::generate_liquid_grid(DensityGrid* terrain_grid) {
    Ref<DensityGrid> liquid_grid;
    liquid_grid.instantiate();
    
    if (!terrain_grid) return liquid_grid;

    int gsx = terrain_grid->get_grid_size_x();
    int gsy = terrain_grid->get_grid_size_y();
    int gsz = terrain_grid->get_grid_size_z();
    
    liquid_grid->initialize_grid(gsx, gsy, gsz, 0.0f);

    Ref<FastNoiseLite> distribution_noise;
    distribution_noise.instantiate();
    distribution_noise->set_seed(noise_seed + 999);
    distribution_noise->set_frequency(0.1f); 
    distribution_noise->set_noise_type(FastNoiseLite::TYPE_PERLIN);

    float terrain_thresh = terrain_grid->get_surface_threshold();
    
    // Helper: Is this voxel "Solid" enough to hold water?
    // We bind these variables by value for lambda use
    auto is_support = [&](Vector3i p) {
        if (p.x < 0 || p.x >= gsx || p.y < 0 || p.y >= gsy || p.z < 0 || p.z >= gsz) return false;
        if (terrain_grid->get_cell(p, 1.0f) >= terrain_thresh) return true; 
        if (liquid_grid->get_cell(p, 0.0f) > 0.5f) return true; 
        return false;
    };
    
    auto is_blocked = [&](Vector3i p) {
        if (p.x < 0 || p.x >= gsx || p.y < 0 || p.y >= gsy || p.z < 0 || p.z >= gsz) return true; 
        if (terrain_grid->get_cell(p, 1.0f) >= terrain_thresh) return true;
        if (liquid_grid->get_cell(p, 0.0f) > 0.5f) return true;
        return false;
    };

    for (int y = 1; y < gsy - 1; ++y) {
        std::unordered_set<int64_t> layer_visited;
        
        for (int z = 0; z < gsz; ++z) {
            for (int x = 0; x < gsx; ++x) {
                Vector3i pos(x, y, z);
                int64_t idx = liquid_grid->get_index(pos);
                
                if (layer_visited.count(idx)) continue; 
                
                if (is_blocked(pos)) continue; 
                if (!is_support(Vector3i(x, y-1, z))) continue; 
                
                if (distribution_noise->get_noise_2d((float)x, (float)z) < sparsity_cutoff) continue;

                std::vector<Vector3i> candidates;
                std::queue<Vector3i> q;
                std::unordered_set<int64_t> current_flood_set; 

                bool leaked = false;
                
                q.push(pos);
                candidates.push_back(pos);
                current_flood_set.insert(idx);
                layer_visited.insert(idx);

                while (!q.empty()) {
                    Vector3i curr = q.front();
                    q.pop();

                    if (candidates.size() > max_basin_size) {
                        leaked = true; 
                        break;
                    }

                    Vector3i nbs[4] = {
                        curr + Vector3i(1,0,0), curr + Vector3i(-1,0,0),
                        curr + Vector3i(0,0,1), curr + Vector3i(0,0,-1)
                    };

                    for (Vector3i nb : nbs) {
                        int64_t nb_idx = liquid_grid->get_index(nb);
                        
                        if (nb.x < 0 || nb.x >= gsx - 1 || nb.z < 0 || nb.z >= gsz - 1 || nb.y < 0 || nb.y >= gsy - 1) {
                            leaked = true;
                            break;
                        }

                        if (is_blocked(nb)) continue;

                        if (!is_support(nb + Vector3i(0, -1, 0))) {
                            leaked = true;
                            break;
                        }

                        if (current_flood_set.find(nb_idx) == current_flood_set.end()) {
                            current_flood_set.insert(nb_idx);
                            layer_visited.insert(nb_idx); 
                            q.push(nb);
                            candidates.push_back(nb);
                        }
                    }
                    if (leaked) break;
                }

                if (!leaked) {
                    for (const Vector3i &p : candidates) {
                        liquid_grid->set_cell(p, water_height_density);

                        Vector3i cardinal_dirs[4] = {
                            Vector3i(1, 0, 0), Vector3i(-1, 0, 0),
                            Vector3i(0, 0, 1), Vector3i(0, 0, -1)
                        };

                        for (const Vector3i &dir : cardinal_dirs) {
                            Vector3i nb = p + dir;
                            if (nb.x < 0 || nb.x >= gsx || nb.z < 0 || nb.z >= gsz) continue;

                            if (terrain_grid->get_cell(nb, 1.0f) >= terrain_thresh) {
                                liquid_grid->set_cell(nb, water_height_density);
                            }
                        }
                    }
                }
            }
        }
    }

    liquid_grid->set_surface_threshold(0.5f);
    return liquid_grid;
}

} // namespace godot