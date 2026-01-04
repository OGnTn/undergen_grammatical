#ifndef CHUNK_MANAGER_H
#define CHUNK_MANAGER_H

#include <godot_cpp/classes/node.hpp>
#include <godot_cpp/variant/typed_array.hpp>
#include "mc_chunk.h"
#include "level_density_grid.h"

namespace godot {

class ChunkManager : public Node {
    GDCLASS(ChunkManager, Node);

protected:
    static void _bind_methods();

public:
    ChunkManager();
    ~ChunkManager();

    // This function takes all necessary data and performs the loop in C++
    void init_chunks_batch(
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
    );
};

} // namespace godot

#endif // CHUNK_MANAGER_H