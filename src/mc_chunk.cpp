// src/mc_chunk.cpp
#include "mc_chunk.h"

#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/variant/utility_functions.hpp>
#include <godot_cpp/core/math.hpp>
#include <godot_cpp/classes/node.hpp> // For notifications like NOTIFICATION_READY
#include <string> // For std::string used in cache key
#include <map>    // For std::map (alternative to Dictionary for vertex cache)
#include <godot_cpp/classes/static_body3d.hpp> // Added for collision body
#include <godot_cpp/classes/array_mesh.hpp> // Added for ArrayMesh type
#include <godot_cpp/classes/occluder_instance3d.hpp> // For occluder node
#include <godot_cpp/classes/array_occluder3d.hpp>   // For occluder shape
#include <godot_cpp/classes/rd_shader_spirv.hpp>
#include <godot_cpp/classes/rd_uniform.hpp>

namespace godot {

MCChunk::MCChunk() {
    // Ensure mesh is set early, even if empty initially
    Ref<ArrayMesh> new_mesh;
    new_mesh.instantiate();
    set_mesh(new_mesh);
}

MCChunk::~MCChunk() {
    // Destructor logic if needed
    _clear_collision(); // Ensure collision is cleaned up on destruction
    _clear_occluder();  // Ensure occluder is cleaned up on destruction
}

void MCChunk::_notification(int p_what) {
    // Optional: If you need behavior like _ready, you can use notifications.
    // For example, generating the mesh when the node enters the tree.
    // if (p_what == NOTIFICATION_READY) {
    //     generate_mesh_from_density_grid();
    // }
}

void MCChunk::_bind_methods() {
    // --- Main Mesh Generation Function ---
    ClassDB::bind_method(D_METHOD("generate_mesh_from_density_grid"), &MCChunk::generate_mesh_from_density_grid);

    // --- Property Getters/Setters ---
    ClassDB::bind_method(D_METHOD("set_chunk_size", "size"), &MCChunk::set_chunk_size);
    ClassDB::bind_method(D_METHOD("get_chunk_size"), &MCChunk::get_chunk_size);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "chunk_size"), "set_chunk_size", "get_chunk_size");

    ClassDB::bind_method(D_METHOD("set_voxel_size", "size"), &MCChunk::set_voxel_size);
    ClassDB::bind_method(D_METHOD("get_voxel_size"), &MCChunk::get_voxel_size);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "voxel_size"), "set_voxel_size", "get_voxel_size");

    ClassDB::bind_method(D_METHOD("set_chunk_grid_offset", "offset"), &MCChunk::set_chunk_grid_offset);
    ClassDB::bind_method(D_METHOD("get_chunk_grid_offset"), &MCChunk::get_chunk_grid_offset);
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3I, "chunk_grid_offset"), "set_chunk_grid_offset", "get_chunk_grid_offset");

    ClassDB::bind_method(D_METHOD("set_terrain_material", "material"), &MCChunk::set_terrain_material);
    ClassDB::bind_method(D_METHOD("get_terrain_material"), &MCChunk::get_terrain_material);
    ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "terrain_material", PROPERTY_HINT_RESOURCE_TYPE, "Material"), "set_terrain_material", "get_terrain_material");

    ClassDB::bind_method(D_METHOD("set_density_grid", "grid"), &MCChunk::set_density_grid);
    ClassDB::bind_method(D_METHOD("get_density_grid"), &MCChunk::get_density_grid);
    // Expose the DensityGrid as a property so it can be assigned in the editor
    ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "density_grid", PROPERTY_HINT_RESOURCE_TYPE, "DensityGrid"), "set_density_grid", "get_density_grid");

    // Bind new collision property
    ClassDB::bind_method(D_METHOD("set_generate_collision", "enable"), &MCChunk::set_generate_collision);
    ClassDB::bind_method(D_METHOD("get_generate_collision"), &MCChunk::get_generate_collision);
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "generate_collision"), "set_generate_collision", "get_generate_collision");

    // Bind new occluder property
    ClassDB::bind_method(D_METHOD("set_generate_occluder", "enable"), &MCChunk::set_generate_occluder);
    ClassDB::bind_method(D_METHOD("get_generate_occluder"), &MCChunk::get_generate_occluder);
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "generate_occluder"), "set_generate_occluder", "get_generate_occluder");

    ClassDB::bind_method(D_METHOD("set_compute_shader", "shader"), &MCChunk::set_compute_shader);
    ClassDB::bind_method(D_METHOD("get_compute_shader"), &MCChunk::get_compute_shader);
    ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "compute_shader", PROPERTY_HINT_RESOURCE_TYPE, "RDShaderFile"), "set_compute_shader", "get_compute_shader");
}

void MCChunk::generate_mesh_from_density_grid() {
    _clear_collision(); // Clear any existing collision shapes
    _clear_occluder();  // Clear any existing occluder nodes

    if (!density_grid.is_valid()) {
        UtilityFunctions::printerr("MCChunk.generate_mesh_from_density_grid: DensityGrid is null! Cannot generate mesh for chunk at offset ", chunk_grid_offset, ".");
        Ref<ArrayMesh> current_mesh = get_mesh();
        if (current_mesh.is_valid()) {
            current_mesh->clear_surfaces(); // Clear any previous mesh
        }
        return;
    }
    if (chunk_size <= 0 || voxel_size <= 0) {
        UtilityFunctions::printerr("MCChunk.generate_mesh_from_density_grid: Invalid ChunkSize (", chunk_size, ") or VoxelSize (", voxel_size, ").");
         Ref<ArrayMesh> current_mesh = get_mesh();
        if (current_mesh.is_valid()) {
            current_mesh->clear_surfaces();
        }
        return;
    }

    // Try Compute Shader first
    if (compute_shader.is_valid()) {
        _generate_mesh_with_compute();
        // _generate_mesh_with_compute handles mesh assignment, collision, and occluder internally
        return;
    }

    //UtilityFunctions::print("Generating mesh for chunk at offset ", chunk_grid_offset);

    // 1. Run Marching Cubes
    Dictionary mc_data = _march_cubes();

    // Check if Dictionary and its keys are valid
    if (mc_data.is_empty() || !mc_data.has("vertices") || !mc_data.has("triangles")) {
        UtilityFunctions::print("MCChunk: _march_cubes returned invalid data for chunk ", chunk_grid_offset);
        Ref<ArrayMesh> current_mesh = get_mesh();
        if (current_mesh.is_valid()) {
            current_mesh->clear_surfaces();
        }
        return;
    }

    PackedVector3Array vertices = mc_data["vertices"];
    PackedInt32Array indices = mc_data["triangles"]; // Should be Int32 for indices

    if (vertices.is_empty() || indices.is_empty()) {
        // UtilityFunctions::print("No mesh data generated for chunk ", chunk_grid_offset);
        Ref<ArrayMesh> current_mesh = get_mesh();
        if (current_mesh.is_valid()) {
            current_mesh->clear_surfaces();
        }
        return;
    }

    // 2. Prepare data for ArrayMesh
    PackedVector3Array normals;
    PackedVector2Array uvs;

    // Simple Planar UV mapping (World XY plane) - Adjust scale and projection as needed
    uvs.resize(vertices.size());
    float uv_scale = voxel_size;
    if (uv_scale <= 0.0) uv_scale = 1.0; // Prevent division by zero
    for (int i = 0; i < vertices.size(); ++i) {
        // Use local vertex positions relative to chunk origin for UVs
        uvs[i] = Vector2(vertices[i].x / uv_scale, vertices[i].y / uv_scale);
    }

    // Calculate Normals (simple method: average face normals)
    normals.resize(vertices.size());
    // Initialize normals to zero vector (important!)
    for (int i = 0; i < normals.size(); ++i) {
        normals[i] = Vector3(0, 0, 0);
    }

    for (int i = 0; i < indices.size(); i += 3) {
        // Add bounds checks for safety
        int i1 = indices[i];
        int i2 = indices[i + 1];
        int i3 = indices[i + 2];

        if (i1 < 0 || i1 >= vertices.size() || i2 < 0 || i2 >= vertices.size() || i3 < 0 || i3 >= vertices.size()) {
             UtilityFunctions::printerr("MCChunk: Invalid index found during normal calculation: ", i1, ", ", i2, ", ", i3);
             continue; // Skip this triangle
        }


        Vector3 v1 = vertices[i1];
        Vector3 v2 = vertices[i2];
        Vector3 v3 = vertices[i3];
        Vector3 face_normal = (v2 - v1).cross(v3 - v1);
        // No need to normalize face_normal here, magnitude contributes to weighting

        // Add face normal to each vertex normal
        normals[i1] += face_normal;
        normals[i2] += face_normal;
        normals[i3] += face_normal;
    }

    // Normalize the vertex normals
    for (int i = 0; i < normals.size(); ++i) {
        normals[i] = normals[i].normalized(); // Normalize the sum
    }


    // 3. Create the Mesh Surface using ArrayMesh
    Array arrays; // Use TypedArray<Array>
    arrays.resize(Mesh::ARRAY_MAX);
    arrays[Mesh::ARRAY_VERTEX] = vertices;
    for(int i = 1; i < Mesh::ARRAY_MAX; i++) {
        arrays[i] = Variant();
    }
    arrays[Mesh::ARRAY_INDEX] = indices;
    arrays[Mesh::ARRAY_NORMAL] = normals;
    //arrays[Mesh::ARRAY_TEX_UV] = uvs;
    //arrays[Mesh::ARRAY_BONES] = PackedFloat32Array(); // Optional: Bone indices if using skinning
    //arrays[Mesh::ARRAY_WEIGHTS] = PackedFloat32Array();

    // Get the mesh resource (should be an ArrayMesh)
    Ref<ArrayMesh> array_mesh = get_mesh();
    //Ref<ArrayMesh> new_mesh;
    //new_mesh.instantiate();
    if (!array_mesh.is_valid()) {
         UtilityFunctions::printerr("MCChunk: Mesh is not a valid ArrayMesh!");
         array_mesh.instantiate(); // Create a new one if it got lost somehow
         set_mesh(array_mesh);
         if(!array_mesh.is_valid()) return; // Still failed
    }

    // Clear previous surfaces and add the new one
    array_mesh->clear_surfaces();
    array_mesh->add_surface_from_arrays(Mesh::PRIMITIVE_TRIANGLES, arrays);

    // 4. Apply Material
    if (terrain_material.is_valid()) {
        // Use set_surface_override_material from MeshInstance3D
        set_surface_override_material(0, terrain_material);
    } else {
        // Warning instead of error, as default material might be okay
        // UtilityFunctions::printerr("MCChunk: TerrainMaterial is not set for chunk ", chunk_grid_offset);
    }

    // 5. Generate Collision!
    if (generate_collision) {
        _generate_collision(vertices, indices);
    }

    // 6. Generate Occluder
    if (generate_occluder) {
        _generate_occluder(vertices, indices);
    }

    //UtilityFunctions::print("Mesh generation complete for chunk at offset ", chunk_grid_offset, ". Verts: ", vertices.size(), ", Tris: ", indices.size() / 3);
}


Dictionary MCChunk::_march_cubes() {
    PackedVector3Array output_vertices;
    PackedInt32Array output_triangles; // Use Int32 for indices

    Dictionary result; // Dictionary to return

    if (!density_grid.is_valid()) {
        UtilityFunctions::printerr("MCChunk._march_cubes: DensityGrid is null!");
        result["vertices"] = output_vertices;
        result["triangles"] = output_triangles;
        return result;
    }

    float surface = density_grid->get_surface_threshold();

    // Vertex sharing optimization: Map edge identifier to vertex index
    // Using Godot Dictionary: Key = String(Vector3i(voxel_base) + "_" + edge_index), Value = int (index in output_vertices)
    Dictionary vertex_cache;

    // Define edge to corner pairs (already in McTables.h)

    // Define the 8 corner offsets relative to the voxel's minimum corner
    const Vector3i corner_offsets[8] = {
        Vector3i(0, 0, 0), Vector3i(1, 0, 0), Vector3i(1, 1, 0), Vector3i(0, 1, 0),
        Vector3i(0, 0, 1), Vector3i(1, 0, 1), Vector3i(1, 1, 1), Vector3i(0, 1, 1)
    };


    // Iterate through each voxel cell within this chunk
    for (int z = 0; z < chunk_size; ++z) {
        for (int y = 0; y < chunk_size; ++y) {
            for (int x = 0; x < chunk_size; ++x) {
                Vector3i local_pos(x, y, z);
                Vector3i world_pos_base = chunk_grid_offset + local_pos; // Min corner of the voxel in world space

                // Sample density values at the 8 corners of the voxel
                float corner_values[8];
                Vector3i world_corner_pos[8]; // World positions of corners

                for (int i = 0; i < 8; ++i) {
                    world_corner_pos[i] = world_pos_base + corner_offsets[i];
                    // Default to solid if outside grid
                    // Use LevelDensityGrid constant here if available, otherwise define locally or get from density_grid if possible
                    corner_values[i] = density_grid->get_cell(world_corner_pos[i], 1.0f);
                }

                // Determine the index based on which corners are "inside"
                int cube_index = 0;
                for (int i = 0; i < 8; ++i) {
                    if (corner_values[i] > surface) { // Compare with surface threshold
                        cube_index |= (1 << i);
                    }
                }

                // Look up the triangles for this cube configuration
                // Access the static data from McTables
                const int* tri_table_row = McTables::TRI_TABLE[cube_index]; // Get pointer to the row

                if (tri_table_row[0] == -1) continue; // Empty voxel

                // Define the 8 corners of the voxel in local chunk space (relative to chunk origin)
                Vector3 corner_locations_local[8];
                for (int i = 0; i < 8; ++i) {
                    // Calculate local corner position based on local_pos and offset, scaled by voxel_size
                    corner_locations_local[i] = Vector3(local_pos + corner_offsets[i]) * voxel_size;
                }

                // Generate vertices and triangles for this voxel
                int triangle_indices[3]; // Store 3 indices per triangle

                // Iterate through the tri_table_row for this cube_index
                for (int i = 0; tri_table_row[i] != -1; i += 3) { // Process 3 edge indices at a time until -1
                    for (int j = 0; j < 3; ++j) { // For each of the 3 vertices in the triangle
                        int edge_index = tri_table_row[i + j];

                        // Create a unique key for the vertex cache based on the edge's world position and index
                        // An edge is uniquely defined by its starting corner (world_pos_base) and the edge index (0-11).
                        // We need to be careful: multiple voxels share edges. The key needs to be consistent.
                        // A common approach is to use the minimum coordinate corner of the edge as the anchor.
                        int corner_a_idx = McTables::CORNER_INDEX_A_FROM_EDGE[edge_index];
                        int corner_b_idx = McTables::CORNER_INDEX_B_FROM_EDGE[edge_index];
                        Vector3i world_corner_a = world_pos_base + corner_offsets[corner_a_idx];
                        Vector3i world_corner_b = world_pos_base + corner_offsets[corner_b_idx];

                        // Determine the canonical edge representation (e.g., based on min coords)
                        // This is complex. A simpler cache key: voxel base + edge index.
                        // This means vertices on shared edges might be duplicated.
                        // Let's stick to the simpler key for now, matching the GDScript version.
                        String edge_full_key = String::num_int64(world_pos_base.x) + "_" +
                                               String::num_int64(world_pos_base.y) + "_" +
                                               String::num_int64(world_pos_base.z) + "_" +
                                               String::num_int64(edge_index);


                        if (vertex_cache.has(edge_full_key)) {
                            triangle_indices[j] = (int)vertex_cache[edge_full_key]; // Cast Variant to int
                        } else {
                            // Calculate vertex if not cached
                            // Indices corner_a_idx and corner_b_idx are correct relative to the current voxel
                            Vector3 vert_pos = _interpolate_vertex(
                                                    corner_locations_local[corner_a_idx],
                                                    corner_locations_local[corner_b_idx],
                                                    corner_values[corner_a_idx],
                                                    corner_values[corner_b_idx]
                                                );

                            // Add vertex to output array and cache its index
                            int new_vertex_index = output_vertices.size();
                            output_vertices.append(vert_pos);
                            vertex_cache[edge_full_key] = new_vertex_index; // Store index in cache
                            triangle_indices[j] = new_vertex_index;
                        }
                    }

                    // Add the triangle indices (ensuring correct winding order - MC tables are usually CCW)
                    output_triangles.append(triangle_indices[0]);
                    output_triangles.append(triangle_indices[1]);
                    output_triangles.append(triangle_indices[2]);
                } // End loop through triangle edges for this voxel
            } // End x loop
        } // End y loop
    } // End z loop


    result["vertices"] = output_vertices;
    result["triangles"] = output_triangles;
    return result;
}


Vector3 MCChunk::_interpolate_vertex(const Vector3 &p1, const Vector3 &p2, float val1, float val2) {
     if (!density_grid.is_valid()) {
        // Fallback if density grid is missing
        return (p1 + p2) * 0.5f;
    }

    float surface = density_grid->get_surface_threshold();

    // Use a small epsilon for float comparison
    if (Math::abs(val1 - val2) < CMP_EPSILON) {
        // Return midpoint if values are too close to avoid division by zero / instability
        return (p1 + p2) * 0.5f;
    }

    // Linear interpolation factor
    // Ensure denominator is not zero before division
    float denominator = val2 - val1;
    // Note: This check might be redundant due to the CMP_EPSILON check above, but added for safety.
    if (Math::abs(denominator) < CMP_EPSILON) {
        return (p1 + p2) * 0.5f;
    }

    float t = (surface - val1) / denominator;

    // Clamp t to avoid extrapolation issues if surface level is outside the range [val1, val2]
    t = Math::clamp(t, 0.0f, 1.0f);

    // Use Godot's lerp function (Vector3 has lerp)
    return p1.lerp(p2, t);
}

void MCChunk::_generate_mesh_with_compute() {
    // 1. Use a LOCAL RenderingDevice.
    // We use a static raw pointer because RenderingDevice is NOT RefCounted.
    static RenderingDevice* local_rd = nullptr;
    
    if (local_rd == nullptr) {
        local_rd = RenderingServer::get_singleton()->create_local_rendering_device();
    }
    
    // Use the local device pointer
    RenderingDevice *rd = local_rd;
    if (!rd) return;

    Ref<RDShaderSPIRV> shader_spirv = compute_shader->get_spirv();
    if (shader_spirv.is_null()) return;

    RID shader_rid = rd->shader_create_from_spirv(shader_spirv);
    if (!shader_rid.is_valid()) return;

    // 1. Prepare Buffers
    
    // A. Density Grid Buffer
    PackedFloat32Array density_data = density_grid->get_world_density_grid();
    PackedByteArray density_bytes = density_data.to_byte_array();
    RID density_buffer = rd->storage_buffer_create(density_bytes.size(), density_bytes);

    // B. Triangulation Table Buffer
    // Flatten the 2D array from McTables
    PackedInt32Array tri_table_flat;
    tri_table_flat.resize(256 * 16);
    for (int i = 0; i < 256; ++i) {
        for (int j = 0; j < 16; ++j) {
            tri_table_flat[i * 16 + j] = McTables::TRI_TABLE[i][j];
        }
    }
    PackedByteArray tri_table_bytes = tri_table_flat.to_byte_array();
    RID tri_table_buffer = rd->storage_buffer_create(tri_table_bytes.size(), tri_table_bytes);

    // C. Output Vertex Buffer
    // Estimate max size: chunk_size^3 * 5 triangles * 3 verts * 3 floats * 4 bytes
    // This is a worst-case allocation.
    int max_verts = chunk_size * chunk_size * chunk_size * 15; 
    int vert_buffer_size = max_verts * 3 * sizeof(float);
    RID vertex_buffer = rd->storage_buffer_create(vert_buffer_size);

    // D. Counter Buffer (Atomic)
    PackedByteArray counter_bytes;
    counter_bytes.resize(4); // uint32
    counter_bytes.encode_u32(0, 0);
    RID counter_buffer = rd->storage_buffer_create(4, counter_bytes);

    // 2. Uniform Sets
    Ref<RDUniform> u_density;
    u_density.instantiate();
    u_density->set_uniform_type(RenderingDevice::UNIFORM_TYPE_STORAGE_BUFFER);
    u_density->set_binding(0);
    u_density->add_id(density_buffer);

    Ref<RDUniform> u_table;
    u_table.instantiate();
    u_table->set_uniform_type(RenderingDevice::UNIFORM_TYPE_STORAGE_BUFFER);
    u_table->set_binding(1);
    u_table->add_id(tri_table_buffer);

    Ref<RDUniform> u_verts;
    u_verts.instantiate();
    u_verts->set_uniform_type(RenderingDevice::UNIFORM_TYPE_STORAGE_BUFFER);
    u_verts->set_binding(2);
    u_verts->add_id(vertex_buffer);

    Ref<RDUniform> u_counter;
    u_counter.instantiate();
    u_counter->set_uniform_type(RenderingDevice::UNIFORM_TYPE_STORAGE_BUFFER);
    u_counter->set_binding(3);
    u_counter->add_id(counter_buffer);

    Array uniforms;
    uniforms.push_back(u_density);
    uniforms.push_back(u_table);
    uniforms.push_back(u_verts);
    uniforms.push_back(u_counter);

    RID uniform_set = rd->uniform_set_create(uniforms, shader_rid, 0);

    // 3. Push Constants
    struct PushConstants {
        int chunk_size_x;
        int chunk_size_y;
        int chunk_size_z;
        int grid_dim_x;
        int grid_dim_y;
        int grid_dim_z;
        int offset_x;
        int offset_y;
        int offset_z;
        float voxel_size;
        float surface_level;
        float padding;
    } push_constants;

    push_constants.chunk_size_x = chunk_size;
    push_constants.chunk_size_y = chunk_size;
    push_constants.chunk_size_z = chunk_size;
    push_constants.grid_dim_x = density_grid->get_grid_size_x();
    push_constants.grid_dim_y = density_grid->get_grid_size_y();
    push_constants.grid_dim_z = density_grid->get_grid_size_z();
    push_constants.offset_x = chunk_grid_offset.x;
    push_constants.offset_y = chunk_grid_offset.y;
    push_constants.offset_z = chunk_grid_offset.z;
    push_constants.voxel_size = voxel_size;
    push_constants.surface_level = density_grid->get_surface_threshold();

    PackedByteArray push_constant_bytes;
    push_constant_bytes.resize(sizeof(PushConstants));
    memcpy(push_constant_bytes.ptrw(), &push_constants, sizeof(PushConstants));

    // 4. Dispatch
    RID pipeline = rd->compute_pipeline_create(shader_rid);
    int64_t compute_list = rd->compute_list_begin();
    rd->compute_list_bind_compute_pipeline(compute_list, pipeline);
    rd->compute_list_bind_uniform_set(compute_list, uniform_set, 0);
    rd->compute_list_set_push_constant(compute_list, push_constant_bytes, sizeof(PushConstants));
    
    // Workgroup size is 4x4x4 = 64.
    int groups_x = (chunk_size + 3) / 4;
    int groups_y = (chunk_size + 3) / 4;
    int groups_z = (chunk_size + 3) / 4;
    rd->compute_list_dispatch(compute_list, groups_x, groups_y, groups_z);
    rd->compute_list_end();

    // 5. Sync and Readback
    rd->submit();
    rd->sync(); // Block until done

    // Read counter
    PackedByteArray counter_output = rd->buffer_get_data(counter_buffer);
    uint32_t vertex_count = counter_output.decode_u32(0);

    if (vertex_count == 0) {
        Ref<ArrayMesh> current_mesh = get_mesh();
        if (current_mesh.is_valid()) current_mesh->clear_surfaces();
        _clear_collision();
        _clear_occluder();
        // Cleanup RIDs
        rd->free_rid(shader_rid);
        rd->free_rid(density_buffer);
        rd->free_rid(tri_table_buffer);
        rd->free_rid(vertex_buffer);
        rd->free_rid(counter_buffer);
        return;
    }

    // Read vertices
    // We only need to read up to vertex_count * 3 floats * 4 bytes
    PackedByteArray vertex_bytes = rd->buffer_get_data(vertex_buffer, 0, vertex_count * 3 * sizeof(float));
    PackedVector3Array raw_vertices;
    raw_vertices.resize(vertex_count);
    const float *float_ptr = reinterpret_cast<const float*>(vertex_bytes.ptr());
    for(uint32_t i=0; i<vertex_count; ++i) {
        raw_vertices[i] = Vector3(float_ptr[i*3], float_ptr[i*3+1], float_ptr[i*3+2]);
    }

    // Cleanup GPU resources
    rd->free_rid(shader_rid);
    rd->free_rid(density_buffer);
    rd->free_rid(tri_table_buffer);
    rd->free_rid(vertex_buffer);
    rd->free_rid(counter_buffer);

    // 6. Post-Process (Weld vertices for smooth normals and indexing)
    PackedVector3Array final_vertices;
    PackedInt32Array final_indices;
    PackedVector3Array final_normals;
    std::map<Vector3, int> vertex_map;

    for (int i = 0; i < raw_vertices.size(); i += 3) {
        Vector3 v1 = raw_vertices[i];
        Vector3 v2 = raw_vertices[i+1];
        Vector3 v3 = raw_vertices[i+2];
        
        // Calculate Face Normal
        Vector3 face_normal = (v2 - v1).cross(v3 - v1); 

        int idx[3];
        Vector3 verts[3] = {v1, v2, v3};

        for(int j=0; j<3; ++j) {
            if(vertex_map.find(verts[j]) == vertex_map.end()) {
                // New unique vertex found
                int new_idx = final_vertices.size();
                final_vertices.append(verts[j]);
                
                // Initialize normal for this NEW vertex
                final_normals.append(Vector3(0,0,0)); 
                
                vertex_map[verts[j]] = new_idx;
                idx[j] = new_idx;
            } else {
                // Existing vertex found
                idx[j] = vertex_map[verts[j]];
            }
            
            final_indices.append(idx[j]);
            
            // Accumulate normal for this existing/new vertex
            final_normals[idx[j]] += face_normal;
        }
    }

    // Normalize normals
    for(int i=0; i<final_normals.size(); ++i) {
        final_normals[i] = final_normals[i].normalized();
    }

    // 7. Create Mesh
    Array arrays;
    arrays.resize(Mesh::ARRAY_MAX);
    arrays[Mesh::ARRAY_VERTEX] = final_vertices;
    arrays[Mesh::ARRAY_INDEX] = final_indices;
    arrays[Mesh::ARRAY_NORMAL] = final_normals;

    // UVs (Planar)
    PackedVector2Array uvs;
    uvs.resize(final_vertices.size());
    float uv_scale = voxel_size > 0 ? voxel_size : 1.0f;
    for(int i=0; i<final_vertices.size(); ++i) {
        uvs[i] = Vector2(final_vertices[i].x / uv_scale, final_vertices[i].y / uv_scale);
    }
    // arrays[Mesh::ARRAY_TEX_UV] = uvs; // Uncomment if needed

    Ref<ArrayMesh> array_mesh = get_mesh();
    if (!array_mesh.is_valid()) {
         array_mesh.instantiate();
         set_mesh(array_mesh);
    }
    array_mesh->clear_surfaces();
    array_mesh->add_surface_from_arrays(Mesh::PRIMITIVE_TRIANGLES, arrays);

    if (terrain_material.is_valid()) {
        set_surface_override_material(0, terrain_material);
    }

    if (generate_collision) {
        _generate_collision(final_vertices, final_indices);
    }
    if (generate_occluder) {
        _generate_occluder(final_vertices, final_indices);
    }
}


// --- Collision Generation Implementation ---
void MCChunk::_generate_collision(const PackedVector3Array &p_vertices, const PackedInt32Array &p_indices) {
    if (p_vertices.is_empty() || p_indices.is_empty()) {
        UtilityFunctions::print("MCChunk._generate_collision: No vertices or indices to generate collision.");
        _clear_collision(); // Ensure no leftover collision
        return;
    }

    // Create a StaticBody3D as a child to hold the collision shape
    // This allows the mesh and its collision to be separate nodes, which is good practice.
    StaticBody3D *static_body = memnew(StaticBody3D);
    add_child(static_body);
    static_body->set_owner(get_owner()); // Important for scene saving/duplication

    // Create the ConcavePolygonShape3D resource
    collision_shape.instantiate();
    if (!collision_shape.is_valid()) {
        UtilityFunctions::printerr("MCChunk._generate_collision: Failed to instantiate ConcavePolygonShape3D.");
        static_body->queue_free(); // Clean up the body
        return;
    }

    // Set the mesh data for the concave polygon shape
    // ConcavePolygonShape3D uses a face list, so we convert the indexed triangle list
    // into a flat list of vertices (triangle by triangle).
    PackedVector3Array collision_face_vertices;
    collision_face_vertices.resize(p_indices.size());

    for (int i = 0; i < p_indices.size(); ++i) {
        int vertex_index = p_indices[i];
        if (vertex_index >= 0 && vertex_index < p_vertices.size()) {
            collision_face_vertices[i] = p_vertices[vertex_index];
        } else {
            UtilityFunctions::printerr("MCChunk._generate_collision: Invalid index in p_indices during collision mesh conversion.");
            // Handle error, maybe clear collision_face_vertices and break
            collision_face_vertices.clear();
            static_body->queue_free();
            collision_shape.unref(); // Release the invalid shape
            return;
        }
    }
    collision_shape->set_faces(collision_face_vertices);


    // Create a CollisionShape3D node and assign the shape
    CollisionShape3D *cs = memnew(CollisionShape3D);
    static_body->add_child(cs);
    cs->set_owner(get_owner()); // Important for scene saving/duplication
    cs->set_shape(collision_shape);

    //UtilityFunctions::print("MCChunk: Generated collision for chunk at offset ", chunk_grid_offset);
}

void MCChunk::_clear_collision() {
    // Iterate through children to find and remove StaticBody3D (and its CollisionShape3D)
    // We iterate backwards to safely remove children while iterating.
    for (int i = get_child_count() - 1; i >= 0; --i) {
        Node *child = get_child(i);
        if (Object::cast_to<StaticBody3D>(child)) {
            // Found a StaticBody3D, remove it and its children (the CollisionShape3D)
            child->queue_free(); // Queue for deletion
            UtilityFunctions::print("MCChunk: Cleared existing collision body for chunk at offset ", chunk_grid_offset);
            // If you only expect one StaticBody3D for collision, you can break here.
            // If multiple might exist, continue the loop.
        }
    }
    collision_shape.unref(); // Also unreference the shape resource itself
}


// --- Occluder Generation Implementation ---
void MCChunk::_generate_occluder(const PackedVector3Array &p_vertices, const PackedInt32Array &p_indices) {
    if (p_vertices.is_empty() || p_indices.is_empty()) {
        UtilityFunctions::print("MCChunk._generate_occluder: No vertices or indices to generate occluder.");
        _clear_occluder(); // Ensure no leftover occluder
        return;
    }

    // Create an OccluderInstance3D as a child
    OccluderInstance3D *occluder_instance = memnew(OccluderInstance3D);
    add_child(occluder_instance);
    occluder_instance->set_owner(get_owner()); // Important for scene saving

    // Create the ArrayOccluder3D resource
    occluder_shape.instantiate();
    if (!occluder_shape.is_valid()) {
        UtilityFunctions::printerr("MCChunk._generate_occluder: Failed to instantiate ArrayOccluder3D.");
        occluder_instance->queue_free(); // Clean up the instance
        return;
    }

    // Set the vertices and indices for the occluder shape.
    // ArrayOccluder3D directly uses vertices and an index array.
    occluder_shape->set_vertices(p_vertices);
    occluder_shape->set_indices(p_indices);

    // Assign the shape to the OccluderInstance3D node
    occluder_instance->set_occluder(occluder_shape);

    //UtilityFunctions::print("MCChunk: Generated occluder for chunk at offset ", chunk_grid_offset);
}

void MCChunk::_clear_occluder() {
    // Iterate through children to find and remove OccluderInstance3D
    // We iterate backwards to safely remove children while iterating.
    for (int i = get_child_count() - 1; i >= 0; --i) {
        Node *child = get_child(i);
        if (Object::cast_to<OccluderInstance3D>(child)) {
            // Found an OccluderInstance3D, remove it
            child->queue_free(); // Queue for deletion
            UtilityFunctions::print("MCChunk: Cleared existing occluder instance for chunk at offset ", chunk_grid_offset);
        }
    }
    occluder_shape.unref(); // Also unreference the shape resource itself
}


// --- Property Getters/Setters Implementation ---

void MCChunk::set_chunk_size(int p_size) { chunk_size = p_size > 0 ? p_size : 1; }
int MCChunk::get_chunk_size() const { return chunk_size; }

void MCChunk::set_voxel_size(float p_size) { voxel_size = p_size > 0 ? p_size : 0.001f; }
float MCChunk::get_voxel_size() const { return voxel_size; }

void MCChunk::set_chunk_grid_offset(const Vector3i &p_offset) { chunk_grid_offset = p_offset; }
Vector3i MCChunk::get_chunk_grid_offset() const { return chunk_grid_offset; }

void MCChunk::set_terrain_material(const Ref<Material> &p_material) { terrain_material = p_material; }
Ref<Material> MCChunk::get_terrain_material() const { return terrain_material; }

void MCChunk::set_density_grid(const Ref<DensityGrid> &p_grid) { density_grid = p_grid; }
Ref<DensityGrid> MCChunk::get_density_grid() const { return density_grid; }

void MCChunk::set_generate_collision(bool p_generate) {
    generate_collision = p_generate;
    // If you want collision to immediately update when this property changes in editor:
    // generate_mesh_from_density_grid(); // This will regenerate mesh AND collision
}
bool MCChunk::get_generate_collision() const {
    return generate_collision;
}

void MCChunk::set_generate_occluder(bool p_generate) {
    generate_occluder = p_generate;
}
bool MCChunk::get_generate_occluder() const {
    return generate_occluder;
}

void MCChunk::set_compute_shader(const Ref<RDShaderFile> &p_shader) { compute_shader = p_shader; }
Ref<RDShaderFile> MCChunk::get_compute_shader() const { return compute_shader; }

} // namespace godot
