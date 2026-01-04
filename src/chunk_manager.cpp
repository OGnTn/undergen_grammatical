#include "chunk_manager.h"
#include <godot_cpp/classes/engine.hpp>
#include <godot_cpp/classes/scene_tree.hpp>
#include <godot_cpp/classes/time.hpp> // Required for profiling
#include <godot_cpp/variant/utility_functions.hpp>

namespace godot {

ChunkManager::ChunkManager() {}
ChunkManager::~ChunkManager() {}

void ChunkManager::_bind_methods() {
    ClassDB::bind_method(D_METHOD("init_chunks_batch", "parent", "chunks_list", "liquid_chunks_list", "density_grid", "terrain_material", "liquid_material", "compute_shader", "chunk_size", "voxel_size", "count_x", "count_y", "count_z"), &ChunkManager::init_chunks_batch);
}

void ChunkManager::init_chunks_batch(
    Node* parent,
    TypedArray<MCChunk> chunks_list,
    TypedArray<MCChunk> liquid_chunks_list,
    Ref<LevelDensityGrid> density_grid,
    TypedArray<Material> terrain_materials,
    Ref<Material> liquid_material,
    Ref<RDShaderFile> compute_shader,
    int chunk_size,
    float voxel_size,
    int count_x, 
    int count_y, 
    int count_z
) {
    if (!parent || !density_grid.is_valid()) {
        UtilityFunctions::printerr("ChunkManager: Invalid parent or density grid.");
        return;
    }

    // --- PROFILING START ---
    Time* time = Time::get_singleton();
    uint64_t start_time = time->get_ticks_msec();
    
    // Accumulators for profiling (in microseconds)
    uint64_t t_alloc = 0;      // Time spent on memnew and variable creation
    uint64_t t_setup = 0;      // Time spent setting properties (materials, size, etc)
    uint64_t t_add_child = 0;  // Time spent interacting with SceneTree
    uint64_t t_set_owner = 0;  // Time spent setting owner (editor only)
    uint64_t t_array = 0;      // Time spent populating the TypedArrays
    
    // --- PREPARATION ---
    int total_chunks = count_x * count_y * count_z;
    
    // Resize arrays (usually fast)
    chunks_list.resize(total_chunks);
    liquid_chunks_list.resize(total_chunks);

    bool is_editor = Engine::get_singleton()->is_editor_hint();
    Node* scene_root = nullptr;
    if (is_editor && parent->get_tree()) {
        scene_root = parent->get_tree()->get_edited_scene_root();
    }

    // --- MAIN LOOP ---
    for (int x = 0; x < count_x; ++x) {
        for (int y = 0; y < count_y; ++y) {
            for (int z = 0; z < count_z; ++z) {
                
                uint64_t t1 = time->get_ticks_usec();
                
                // 1. ALLOCATION
                int index = x + y * count_x + z * count_x * count_y;
                Vector3i grid_offset(x, y, z);
                Vector3 pos(x * chunk_size * voxel_size, y * chunk_size * voxel_size, z * chunk_size * voxel_size);

                MCChunk* chunk = memnew(MCChunk);
                String t_name = "TerrainChunk_" + String::num_int64(x) + "_" + String::num_int64(y) + "_" + String::num_int64(z);
                
                uint64_t t2 = time->get_ticks_usec();
                t_alloc += (t2 - t1);

                // 2. SETUP PROPERTIES
                chunk->set_name(t_name);
                chunk->set_chunk_size(chunk_size);
                chunk->set_voxel_size(voxel_size);
                chunk->set_chunk_grid_offset(grid_offset * chunk_size);
                chunk->set_materials(terrain_materials);
                chunk->set_density_grid(density_grid);
                chunk->set_generate_collision(true);
                chunk->set_compute_shader(compute_shader);
                chunk->set_position(pos);
                
                uint64_t t3 = time->get_ticks_usec();
                t_setup += (t3 - t2);

                // 3. ADD CHILD (Likely the bottleneck)
                parent->add_child(chunk); // true = force readable name
                
                uint64_t t4 = time->get_ticks_usec();
                t_add_child += (t4 - t3);

                // 4. SET OWNER
                if (is_editor && scene_root) {
                    chunk->set_owner(scene_root);
                }
                
                uint64_t t5 = time->get_ticks_usec();
                t_set_owner += (t5 - t4);

                // 5. ARRAY ASSIGNMENT
                chunks_list[index] = chunk;
                
                uint64_t t6 = time->get_ticks_usec();
                t_array += (t6 - t5);

                // --- REPEAT FOR LIQUID (Simplified for brevity, accumulators continue) ---
                
                uint64_t l_t1 = time->get_ticks_usec();
                MCChunk* liquid_chunk = memnew(MCChunk);
                String l_name = "LiquidChunk_" + String::num_int64(x) + "_" + String::num_int64(y) + "_" + String::num_int64(z);
                uint64_t l_t2 = time->get_ticks_usec();
                t_alloc += (l_t2 - l_t1);

                liquid_chunk->set_name(l_name);
                liquid_chunk->set_chunk_size(chunk_size);
                int res_mult = density_grid->get_liquid_resolution_multiplier();
                liquid_chunk->set_voxel_size(voxel_size / (float)res_mult);
                liquid_chunk->set_chunk_grid_offset(grid_offset * chunk_size);
                TypedArray<Material> liquid_mats;
                liquid_mats.append(liquid_material);
                liquid_chunk->set_materials(liquid_mats);
                liquid_chunk->set_generate_collision(false);
                liquid_chunk->set_generate_occluder(false);
                liquid_chunk->set_position(pos);
                uint64_t l_t3 = time->get_ticks_usec();
                t_setup += (l_t3 - l_t2);

                parent->add_child(liquid_chunk);
                uint64_t l_t4 = time->get_ticks_usec();
                t_add_child += (l_t4 - l_t3);

                if (is_editor && scene_root) liquid_chunk->set_owner(scene_root);
                uint64_t l_t5 = time->get_ticks_usec();
                t_set_owner += (l_t5 - l_t4);

                liquid_chunks_list[index] = liquid_chunk;
                uint64_t l_t6 = time->get_ticks_usec();
                t_array += (l_t6 - l_t5);
            }
        }
    }
    
    uint64_t end_time = time->get_ticks_msec();

    // --- PRINT REPORT ---
    UtilityFunctions::print("=== Chunk Initialization Profiling ===");
    UtilityFunctions::print("Total Time: ", end_time - start_time, " ms");
    UtilityFunctions::print("Total Chunks Created: ", total_chunks * 2); // Terrain + Liquid
    UtilityFunctions::print("-------------------------------------");
    UtilityFunctions::print("Allocation (memnew):  ", t_alloc / 1000, " ms");
    UtilityFunctions::print("Property Setup:       ", t_setup / 1000, " ms");
    UtilityFunctions::print("SceneTree (add_child):", t_add_child / 1000, " ms");
    UtilityFunctions::print("Set Owner (Editor):   ", t_set_owner / 1000, " ms");
    UtilityFunctions::print("Array Assignment:     ", t_array / 1000, " ms");
    UtilityFunctions::print("======================================");
}

} // namespace godot