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
#include <godot_cpp/classes/time.hpp>

namespace godot {

struct SurfaceBuilder {
    PackedVector3Array vertices;
    PackedVector3Array normals;
    PackedInt32Array indices;
    // Maps "Master Index" -> "Surface Index" to weld vertices within the surface
    std::map<int, int> index_cache; 
};

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

    ClassDB::bind_method(D_METHOD("set_materials", "materials"), &MCChunk::set_materials);
    ClassDB::bind_method(D_METHOD("get_materials"), &MCChunk::get_materials);
    ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "materials", PROPERTY_HINT_TYPE_STRING, "4/17:Material"), "set_materials", "get_materials");

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
    _clear_collision();
    _clear_occluder();
    
    if (!density_grid.is_valid()) {
        Ref<ArrayMesh> current_mesh = get_mesh();
        if (current_mesh.is_valid()) current_mesh->clear_surfaces();
        return;
    }

    Time* time = Time::get_singleton();
    uint64_t start_total = time->get_ticks_usec();
    uint64_t t1 = start_total;

    // --- STEP 1: Generate Unified Geometry (The "Master" Mesh) ---
    // We collect ALL triangles first, regardless of material, so we can smooth the normals.
    
    PackedVector3Array master_vertices;
    PackedInt32Array master_indices;
    std::vector<int> master_triangle_materials; // Stores the material ID for each triangle
    
    // Vertex Cache for the master mesh (Edge Key -> Master Index)
    std::map<String, int> master_vertex_cache;

    float surface = density_grid->get_surface_threshold();
    const Vector3i corner_offsets[8] = {
        Vector3i(0, 0, 0), Vector3i(1, 0, 0), Vector3i(1, 1, 0), Vector3i(0, 1, 0),
        Vector3i(0, 0, 1), Vector3i(1, 0, 1), Vector3i(1, 1, 1), Vector3i(0, 1, 1)
    };

    // Run Marching Cubes
    for (int z = 0; z < chunk_size; ++z) {
        for (int y = 0; y < chunk_size; ++y) {
            for (int x = 0; x < chunk_size; ++x) {
                Vector3i local_pos(x, y, z);
                Vector3i world_pos_base = chunk_grid_offset + local_pos;

                // 1. Get Material ID for this voxel
                int mat_idx = _get_voxel_material_id(local_pos);

                // 2. Sample Corners
                float corner_values[8];
                for (int i = 0; i < 8; ++i) {
                    corner_values[i] = density_grid->get_cell(world_pos_base + corner_offsets[i], 1.0f);
                }

                // 3. Determine Cube Index
                int cube_index = 0;
                for (int i = 0; i < 8; ++i) {
                    if (corner_values[i] > surface) cube_index |= (1 << i);
                }

                const int* tri_table_row = McTables::TRI_TABLE[cube_index];
                if (tri_table_row[0] == -1) continue;

                // 4. Precompute Local Corner Positions
                Vector3 corner_locs[8];
                for (int i = 0; i < 8; ++i) {
                    corner_locs[i] = Vector3(local_pos + corner_offsets[i]) * voxel_size;
                }

                // 5. Generate Triangles
                for (int i = 0; tri_table_row[i] != -1; i += 3) {
                    int tri_indices[3];

                    // Process 3 vertices for the triangle
                    for (int j = 0; j < 3; ++j) {
                        int edge_index = tri_table_row[i + j];
                        
                        // Create unique key for this specific edge in world space
                        String edge_key = String::num_int64(world_pos_base.x) + "_" +
                                          String::num_int64(world_pos_base.y) + "_" +
                                          String::num_int64(world_pos_base.z) + "_" +
                                          String::num_int64(edge_index);

                        // Check Cache
                        auto cache_it = master_vertex_cache.find(edge_key);
                        if (cache_it != master_vertex_cache.end()) {
                            tri_indices[j] = cache_it->second;
                        } else {
                            // Interpolate new vertex
                            int c_a = McTables::CORNER_INDEX_A_FROM_EDGE[edge_index];
                            int c_b = McTables::CORNER_INDEX_B_FROM_EDGE[edge_index];
                            
                            Vector3 vert_pos = _interpolate_vertex(
                                corner_locs[c_a], corner_locs[c_b],
                                corner_values[c_a], corner_values[c_b]
                            );

                            int new_idx = master_vertices.size();
                            master_vertices.append(vert_pos);
                            master_vertex_cache[edge_key] = new_idx;
                            tri_indices[j] = new_idx;
                        }
                    }

                    // Store Triangle Data
                    master_indices.append(tri_indices[0]);
                    master_indices.append(tri_indices[1]);
                    master_indices.append(tri_indices[2]);
                    
                    // Store the material this triangle belongs to
                    master_triangle_materials.push_back(mat_idx);
                }
            }
        }
    }
    
    uint64_t t2 = time->get_ticks_usec();
    uint64_t dt_march = t2 - t1;

    if (master_vertices.is_empty()) return;

    // --- STEP 2: Calculate Smooth Normals on Unified Mesh ---
    PackedVector3Array master_normals;
    master_normals.resize(master_vertices.size());
    // Zero out normals
    for(int i=0; i<master_normals.size(); ++i) master_normals[i] = Vector3(0,0,0);

    // Accumulate face normals
    for (int i = 0; i < master_indices.size(); i += 3) {
        int i1 = master_indices[i];
        int i2 = master_indices[i+1];
        int i3 = master_indices[i+2];

        Vector3 v1 = master_vertices[i1];
        Vector3 v2 = master_vertices[i2];
        Vector3 v3 = master_vertices[i3];

        Vector3 face_normal = (v2 - v1).cross(v3 - v1);
        
        master_normals[i1] += face_normal;
        master_normals[i2] += face_normal;
        master_normals[i3] += face_normal;
    }

    // Normalize
    for(int i=0; i<master_normals.size(); ++i) {
        master_normals[i] = master_normals[i].normalized();
    }
    
    uint64_t t3 = time->get_ticks_usec();
    uint64_t dt_normals = t3 - t2;

    // --- STEP 3: Split into Material Surfaces ---
    std::map<int, SurfaceBuilder> surfaces;

    int num_triangles = master_triangle_materials.size();
    
    for(int t = 0; t < num_triangles; ++t) {
        int mat_id = master_triangle_materials[t];
        SurfaceBuilder &surf = surfaces[mat_id];

        // Process the 3 indices of this triangle
        for(int k = 0; k < 3; ++k) {
            int master_idx = master_indices[t * 3 + k];
            
            // Local Welding: Check if this master vertex is already in THIS surface
            auto it = surf.index_cache.find(master_idx);
            if (it != surf.index_cache.end()) {
                // Reuse existing vertex in this surface
                surf.indices.append(it->second);
            } else {
                // Add new vertex to this surface, copying data from Master
                int new_surf_idx = surf.vertices.size();
                
                surf.vertices.append(master_vertices[master_idx]);
                surf.normals.append(master_normals[master_idx]); // COPY THE SMOOTH NORMAL
                
                surf.indices.append(new_surf_idx);
                surf.index_cache[master_idx] = new_surf_idx;
            }
        }
    }
    
    uint64_t t4 = time->get_ticks_usec();
    uint64_t dt_split = t4 - t3;

    // --- STEP 4: Build Godot Mesh ---
    Ref<ArrayMesh> array_mesh = get_mesh();
    if (!array_mesh.is_valid()) {
         array_mesh.instantiate();
         set_mesh(array_mesh);
    }
    array_mesh->clear_surfaces();

    // Prepare collision data (collects all verts again)
    PackedVector3Array total_collision_vertices;
    PackedInt32Array total_collision_indices;
    int col_offset = 0;

    // Iterate over our split surfaces
    for (auto const& [mat_id, builder] : surfaces) {
        if (builder.vertices.is_empty()) continue;

        Array arrays;
        arrays.resize(Mesh::ARRAY_MAX);
        arrays[Mesh::ARRAY_VERTEX] = builder.vertices;
        arrays[Mesh::ARRAY_NORMAL] = builder.normals;
        arrays[Mesh::ARRAY_INDEX] = builder.indices;

        // Generate Planar UVs
        PackedVector2Array uvs;
        uvs.resize(builder.vertices.size());
        float uv_scale = voxel_size > 0 ? voxel_size : 1.0f;
        for(int i=0; i < builder.vertices.size(); ++i) {
            uvs[i] = Vector2(builder.vertices[i].x / uv_scale, builder.vertices[i].y / uv_scale);
        }
        arrays[Mesh::ARRAY_TEX_UV] = uvs;

        // Create Surface
        array_mesh->add_surface_from_arrays(Mesh::PRIMITIVE_TRIANGLES, arrays);

        // Assign Material
        if (mat_id >= 0 && mat_id < materials.size()) {
            Ref<Material> mat = materials[mat_id];
            if(mat.is_valid()) {
                set_surface_override_material(array_mesh->get_surface_count() - 1, mat);
            }
        }

        // Collect Collision Data
        if (generate_collision || generate_occluder) {
            total_collision_vertices.append_array(builder.vertices);
            for(int idx : builder.indices) {
                total_collision_indices.append(idx + col_offset);
            }
            col_offset += builder.vertices.size();
        }
    }
    
    uint64_t t5 = time->get_ticks_usec();
    uint64_t dt_mesh = t5 - t4;

    // --- STEP 5: Generate Physics/Occlusion ---
    if (generate_collision) _generate_collision(total_collision_vertices, total_collision_indices);
    if (generate_occluder) _generate_occluder(total_collision_vertices, total_collision_indices);
    
    uint64_t t6 = time->get_ticks_usec();
    uint64_t dt_collision = t6 - t5;
    
    double total_ms = (t6 - start_total) / 1000.0;
    
    // Only print if it took significant time to avoid spam
    //if (total_ms > 1.0) {
    //    UtilityFunctions::print("MCChunk::generate_mesh [", get_name(), "] Total: ", total_ms, " ms");
    //    UtilityFunctions::print("  March: ", dt_march/1000.0, " ms, Normals: ", dt_normals/1000.0, " ms, Split: ", dt_split/1000.0, " ms, Mesh: ", dt_mesh/1000.0, " ms, Col: ", dt_collision/1000.0, " ms");
    //}
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

int MCChunk::_get_voxel_material_id(const Vector3i &local_pos) {
    if (!density_grid.is_valid()) return 0;

    // Convert local chunk position to global grid position
    // chunk_grid_offset is already in voxel coordinates (e.g., 0, 16, 32) based on your implementation
    Vector3i global_pos = chunk_grid_offset + local_pos;

    // Query the DensityGrid
    // Since get_material_id handles bounds checking internally, this is safe
    return density_grid->get_material_id(global_pos);
}

Dictionary MCChunk::_march_cubes_multi_mat() {
    // Map material index -> SurfaceData
    std::map<int, SurfaceData> surfaces;
    
    if (!density_grid.is_valid()) return Dictionary();
    
    float surface = density_grid->get_surface_threshold();

    // Corner offsets (same as before)
    const Vector3i corner_offsets[8] = {
        Vector3i(0, 0, 0), Vector3i(1, 0, 0), Vector3i(1, 1, 0), Vector3i(0, 1, 0),
        Vector3i(0, 0, 1), Vector3i(1, 0, 1), Vector3i(1, 1, 1), Vector3i(0, 1, 1)
    };

    for (int z = 0; z < chunk_size; ++z) {
        for (int y = 0; y < chunk_size; ++y) {
            for (int x = 0; x < chunk_size; ++x) {
                Vector3i local_pos(x, y, z);
                Vector3i world_pos_base = chunk_grid_offset + local_pos;

                // --- Material Logic ---
                // Determine which material this voxel belongs to.
                // Note: In MC, a triangle might span two voxels. 
                // We typically assign the triangle to the voxel containing the generated geometry.
                int mat_idx = _get_voxel_material_id(local_pos);
                
                // Get or create the surface data for this material
                SurfaceData &current_surface = surfaces[mat_idx];

                // --- Standard MC Sampling (Same as your code) ---
                float corner_values[8];
                for (int i = 0; i < 8; ++i) {
                    corner_values[i] = density_grid->get_cell(world_pos_base + corner_offsets[i], 1.0f);
                }

                int cube_index = 0;
                for (int i = 0; i < 8; ++i) {
                    if (corner_values[i] > surface) cube_index |= (1 << i);
                }

                const int* tri_table_row = McTables::TRI_TABLE[cube_index];
                if (tri_table_row[0] == -1) continue;

                Vector3 corner_locations_local[8];
                for (int i = 0; i < 8; ++i) {
                    corner_locations_local[i] = Vector3(local_pos + corner_offsets[i]) * voxel_size;
                }

                // --- Triangle Generation ---
                for (int i = 0; tri_table_row[i] != -1; i += 3) {
                    int triangle_indices[3];
                    
                    for (int j = 0; j < 3; ++j) {
                        int edge_index = tri_table_row[i + j];
                        
                        // NOTE: To support multiple surfaces, each surface must have its own 
                        // vertex cache. We cannot share vertices between different surfaces/draw calls.
                        String edge_full_key = String::num_int64(world_pos_base.x) + "_" +
                                               String::num_int64(world_pos_base.y) + "_" +
                                               String::num_int64(world_pos_base.z) + "_" +
                                               String::num_int64(edge_index);

                        if (current_surface.vertex_cache.has(edge_full_key)) {
                            triangle_indices[j] = (int)current_surface.vertex_cache[edge_full_key];
                        } else {
                            int corner_a_idx = McTables::CORNER_INDEX_A_FROM_EDGE[edge_index];
                            int corner_b_idx = McTables::CORNER_INDEX_B_FROM_EDGE[edge_index];

                            Vector3 vert_pos = _interpolate_vertex(
                                corner_locations_local[corner_a_idx],
                                corner_locations_local[corner_b_idx],
                                corner_values[corner_a_idx],
                                corner_values[corner_b_idx]
                            );

                            int new_vertex_index = current_surface.vertices.size();
                            current_surface.vertices.append(vert_pos);
                            current_surface.vertex_cache[edge_full_key] = new_vertex_index;
                            triangle_indices[j] = new_vertex_index;
                        }
                    }

                    current_surface.indices.append(triangle_indices[0]);
                    current_surface.indices.append(triangle_indices[1]);
                    current_surface.indices.append(triangle_indices[2]);
                }
            }
        }
    }

    // Convert std::map to Godot Dictionary for return
    Dictionary result;
    for (auto const& [mat_idx, data] : surfaces) {
        Dictionary surface_dict;
        surface_dict["vertices"] = data.vertices;
        surface_dict["indices"] = data.indices;
        result[mat_idx] = surface_dict;
    }
    return result;
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

    if (materials.size() > 0) {
        Ref<Material> mat = materials[0];
        if (mat.is_valid()) {
            set_surface_override_material(0, mat);
        }
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

void MCChunk::set_materials(const TypedArray<Material> &p_materials) { materials = p_materials; }
TypedArray<Material> MCChunk::get_materials() const { return materials; }

void MCChunk::set_compute_shader(const Ref<RDShaderFile> &p_shader) { compute_shader = p_shader; }
Ref<RDShaderFile> MCChunk::get_compute_shader() const { return compute_shader; }

} // namespace godot
