#ifndef LIQUID_GENERATOR_H
#define LIQUID_GENERATOR_H

#include "density_grid.h"
#include <godot_cpp/classes/fast_noise_lite.hpp>
#include <godot_cpp/classes/ref.hpp>

namespace godot {

class LiquidGenerator {
public:
    // Configuration
    int max_basin_size = 40;       
    float sparsity_cutoff = 0.2f;  
    float water_height_density = 0.5f; 
    int liquid_resolution_multiplier = 2; 

    // Dependencies
    int noise_seed = 1337; // Used for distribution noise
    
    // Logic
    Ref<DensityGrid> generate_liquid_grid(DensityGrid* terrain_grid);

private:
    float WORLD_SOLID_VALUE = 1.0f; 
    
    // Note: _is_space_free logic was in previous code but seemed unused by generate_liquid_grid directly?
    // Actually generate_liquid_grid implementation in original code did NOT use `_is_space_free`.
    // It used lambda `is_support` and `is_blocked`.
    // `_is_space_free` was likely for Entity placement (not Liquid generation).
    // I will extract `generate_liquid_grid` logic.
};

} // namespace godot

#endif // LIQUID_GENERATOR_H