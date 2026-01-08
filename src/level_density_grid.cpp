#define OGT_VOX_IMPLEMENTATION
#include "level_density_grid.h"

#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/variant/utility_functions.hpp>
#include <godot_cpp/classes/file_access.hpp>
#include <godot_cpp/core/math.hpp> 
#include <godot_cpp/classes/time.hpp>

// Define missing OGT macros if needed or just use default
// The library handles it.

#include <godot_cpp/classes/multi_mesh.hpp>
#include <godot_cpp/variant/transform3d.hpp>

namespace godot {

LevelDensityGrid::LevelDensityGrid() {
    // Initialize Ref<> types
    noise_generator.instantiate(); 
    rng.instantiate();
    rng->randomize(); 

    if (noise_generator.is_valid()) {
        noise_generator->set_noise_type(FastNoiseLite::TYPE_PERLIN);
        noise_generator->set_seed(noise_seed);
        noise_generator->set_frequency(noise_frequency);
    }
}

LevelDensityGrid::~LevelDensityGrid() {
}

void LevelDensityGrid::generate_level_data(const Vector3i &world_grid_dimensions, float voxel_size, int64_t seed) {
    if (world_grid_dimensions.x <= 0 || world_grid_dimensions.y <= 0 || world_grid_dimensions.z <= 0) {
        UtilityFunctions::printerr("LevelDensityGrid.generate_level_data: Invalid World Grid Dimensions ", world_grid_dimensions);
        return;
    }

    UtilityFunctions::print("Generating Level Data for Grid Size: ", world_grid_dimensions, " with seed: ", seed);

    Time* time = Time::get_singleton();
    uint64_t start_total = time->get_ticks_msec();
    
    uint64_t t_init_grid = 0;
    uint64_t t_rooms_paths = 0;
    uint64_t t_smoothing = 0;
    uint64_t t_normals = 0;

    // 1. Initialize & Seed
    uint64_t t1 = time->get_ticks_usec();
    
    set_noise_seed(seed);
    if (!rng.is_valid()) {
        rng.instantiate();
        UtilityFunctions::printerr("LevelDensityGrid: RNG was invalid, re-initialized.");
    }
    rng->set_seed(seed);

    initialize_grid(world_grid_dimensions.x, world_grid_dimensions.y, world_grid_dimensions.z, WORLD_SOLID_VALUE);
    
    uint64_t t2 = time->get_ticks_usec();
    t_init_grid = (t2 - t1);

    // 2. Generate Rooms (Delegated)
    int placed_count = 0;
    
    // Pass config to generator (or generator uses its internal config? We need to sync config)
    // Actually, our Getters/Setters update the component config directly. 
    // EXCEPT: if properties are not synced?
    // In set_room_count, we do `room_gen.room_count = p_count`.
    // So the component state is already current.
    
    // Note: RoomGenerator needs `DensityGrid*`, `RNG*`.
    TypedArray<Dictionary> generated_rooms = room_gen.generate_rooms(this, rng.ptr(), placed_count);
    
    // 3. Connect Rooms (Delegated)
    path_carver.connect_rooms(
        this, 
        rng.ptr(), 
        noise_generator, 
        generated_rooms, 
        voxel_size, 
        calculated_spawn_position, 
        calculated_end_position
    );
    
    // Populate resolved_rooms local cache for 'get_room_type_at' if needed?
    // generated_rooms is TypedArray<Dictionary>. resolved_rooms is std::vector<ResolvedRoom>.
    // To support `get_room_type_at`, we need to populate `resolved_rooms`.
    resolved_rooms.clear();
    for(int i=0; i<generated_rooms.size(); ++i) {
        Dictionary d = generated_rooms[i];
        ResolvedRoom r;
        r.position = d["start"];
        r.size = Vector3i(d["end"]) - r.position;
        r.type = "generic"; // or from dict if we stored it
        resolved_rooms.push_back(r);
    }
    
    uint64_t t3 = time->get_ticks_usec();
    t_rooms_paths = (t3 - t2);

    // 4. Apply Noise (Internal)
    _apply_noise();

    // 5. Optional Smoothing
    // 5. Optional Smoothing
    if (smooth_terrain && !path_carver.use_square_brush) {
        UtilityFunctions::print("Smoothing terrain with strength: ", smoothing_strength);
        low_pass();
    }
    
    uint64_t t4 = time->get_ticks_usec();
    t_smoothing = (t4 - t3);

    // 6. Calculate Surface Normals
    _calculate_surface_normals();
    
    uint64_t t5 = time->get_ticks_usec();
    t_normals = (t5 - t4);

    uint64_t end_total = time->get_ticks_msec();

    UtilityFunctions::print("=== Generation Profiling ===");
    UtilityFunctions::print("Total Time:      ", end_total - start_total, " ms");
    UtilityFunctions::print("---------------------------");
    UtilityFunctions::print("Grid Init:       ", t_init_grid / 1000, " ms");
    UtilityFunctions::print("Rooms & Paths:   ", t_rooms_paths / 1000, " ms");
    UtilityFunctions::print("Smoothing:       ", t_smoothing / 1000, " ms");
    UtilityFunctions::print("Surface Normals: ", t_normals / 1000, " ms");
    UtilityFunctions::print("===========================");
}

void LevelDensityGrid::generate_from_graph(const TypedArray<Dictionary> &node_data, const TypedArray<Dictionary> &edge_data, float voxel_size, int64_t seed, Dictionary settings) {
    // 1. Initialize
    initialize_grid(get_grid_size_x(), get_grid_size_y(), get_grid_size_z(), WORLD_SOLID_VALUE);
    resolved_rooms.clear();
    current_vox_spawns.clear();
    
    // Cache Maps
    current_vox_spawns.clear();
    current_vox_spawn_map = settings.get("vox_spawn_map", Dictionary());
    current_vox_material_map = settings.get("vox_material_map", Dictionary());

    current_carving_zone_id = 0; 
    rng->set_seed(seed);
    
// 2. Parse Nodes
    std::vector<ResolvedRoom> processing_rooms;
    std::map<String, int> id_to_index;

    for(int i=0; i < node_data.size(); ++i) {
        Dictionary node = node_data[i];
        ResolvedRoom r;
        r.id = node.get("id", "room");
        r.type = node.get("type", "generic");
        r.vox_path = node.get("vox_path", "");
        
        Vector3i min_s = node.get("min_size", room_gen.min_room_size);
        Vector3i max_s = node.get("max_size", room_gen.max_room_size);
        r.size = Vector3i(
            rng->randi_range(min_s.x, max_s.x),
            rng->randi_range(min_s.y, max_s.y),
            rng->randi_range(min_s.z, max_s.z)
        );

        Vector3i grid_center = get_grid_dimensions() / 2;
        Vector3i max_pos = get_grid_dimensions() - r.size;
        r.position = Vector3i(
            rng->randi_range(0, MAX(0, max_pos.x)),
            rng->randi_range(0, MAX(0, max_pos.y)),
            rng->randi_range(0, MAX(0, max_pos.z))
        );
        
        // Parse Constraints
        Dictionary constraints = node.get("constraints", Dictionary());
        if (constraints.has("fixed_pos")) {
             Vector3 fix_v = constraints["fixed_pos"];
             r.is_fixed = true;
             r.fixed_position = Vector3(grid_center) + fix_v;
             r.position = Vector3i(r.fixed_position) - (r.size / 2);
        }
        if (constraints.has("relative_to")) {
            r.has_relative_constraint = true;
            r.relative_target_id = (String)constraints["relative_to"];
            r.relative_offset = constraints.get("offset", Vector3(0,0,0));
        }

        processing_rooms.push_back(r);
        id_to_index[r.id] = i;
    }

    // 3. Parse Edges
    std::vector<ResolvedEdge> edges;
    for(int i=0; i < edge_data.size(); ++i) {
        Dictionary edge_dict = edge_data[i];
        String from_id = edge_dict.get("from", "");
        String to_id = edge_dict.get("to", "");
        String type = edge_dict.get("type", "corridor"); 

        if(id_to_index.count(from_id) && id_to_index.count(to_id)) {
            ResolvedEdge e;
            e.from_index = id_to_index[from_id];
            e.to_index = id_to_index[to_id];
            e.type = type;
            edges.push_back(e);
        }
    }
    
    // --- Guide Curve & Volumes Prep ---
    PackedVector3Array guide_points = settings.get("guide_curve", PackedVector3Array());
    TypedArray<Dictionary> macro_volumes = settings.get("macro_volumes", TypedArray<Dictionary>());
    
    // Calculate Node Depth
    std::vector<int> node_depths(processing_rooms.size(), 0);
    int max_depth = 1;

    if (guide_points.size() > 1 && processing_rooms.size() > 0) {
        std::vector<int> q;
        q.push_back(0); 
        std::vector<bool> visited(processing_rooms.size(), false);
        visited[0] = true;
        
        std::vector<std::vector<int>> adj(processing_rooms.size());
        for(const auto& e : edges) {
            adj[e.from_index].push_back(e.to_index);
            adj[e.to_index].push_back(e.from_index);
        }
        
        int head = 0;
        while(head < q.size()) {
            int u = q[head++];
            for(int v : adj[u]) {
                if(!visited[v]) {
                    visited[v] = true;
                    node_depths[v] = node_depths[u] + 1;
                    if(node_depths[v] > max_depth) max_depth = node_depths[v];
                    q.push_back(v);
                }
            }
        }
    }
    
    // 4. Place Rooms
    if (placement_algorithm == ALGO_FORCE_DIRECTED) {
        _place_rooms_force_directed(processing_rooms, edges, node_depths, max_depth, macro_volumes, guide_points);
    } else {
        _place_rooms_bsp(processing_rooms);
    }
    
    this->resolved_rooms = processing_rooms; 

    // 5. Carve Generic Rooms (Pass 1)
    for (const auto& room : processing_rooms) {
        if (!room.vox_path.is_empty()) continue; // Skip Vox rooms for later
        
        // Carve Generic Room
        std::vector<ResolvedRoom> single_room_vec;
        single_room_vec.push_back(room);
        room_gen.create_rooms_from_data(this, single_room_vec, current_carving_zone_id);
    }

    // 6. Carve Connections
    path_carver.create_paths_from_edges(this, rng.ptr(), noise_generator, processing_rooms, edges);

    // 7. Stamp Vox Rooms (Pass 2 - Overwrite Paths)
    for (const auto& room : processing_rooms) {
         if (room.vox_path.is_empty()) continue; // Skip Generic rooms
         
         UtilityFunctions::print("DEBUG: Processing Vox Room [", room.id, "] with path: '", room.vox_path, "'");

         const ogt_vox_scene* scene = _load_vox(room.vox_path);
         if (scene) {
             UtilityFunctions::print("DEBUG: Scene loaded successfully. Instances: ", scene->num_instances);
             _stamp_vox(room, scene);
         } else {
             UtilityFunctions::printerr("DEBUG: Failed to load scene for room [", room.id, "]");
             // Fallback if loading failed
             std::vector<ResolvedRoom> single_room_vec;
             single_room_vec.push_back(room);
             room_gen.create_rooms_from_data(this, single_room_vec, current_carving_zone_id); 
         }
    }

    // 7. Post-Processing
    if(smooth_terrain) low_pass();
    _calculate_surface_normals();

    // 8. Determine Start/End Positions
    for (const auto& room : processing_rooms) {
        String type_lower = room.type.to_lower();
        
        if (type_lower.contains("start")) {
            // Default to center
            calculated_spawn_position = room.center();
            
            // Check for precise VOX spawn
            if (!room.vox_path.is_empty()) {
                // Find a spawn point within this room's bounds
                // Since current_vox_spawns are global grid coords, we just check inclusion (or just match room ID if we tracked it, but spatial check is fine)
                for (const auto& spawn : current_vox_spawns) {
                    Vector3 p = spawn.position;
                    if (p.x >= room.position.x && p.x < room.position.x + room.size.x &&
                        p.y >= room.position.y && p.y < room.position.y + room.size.y &&
                        p.z >= room.position.z && p.z < room.position.z + room.size.z) {
                        calculated_spawn_position = p;
                        break; // Found one
                    }
                }
            }
            UtilityFunctions::print("LevelDensityGrid: Start Position set to ", calculated_spawn_position, " (Room: ", room.id, ")");
        } 
        else if (type_lower.contains("end")) {
            // Default to center
            calculated_end_position = room.center();
             
             // Check for precise VOX spawn
            if (!room.vox_path.is_empty()) {
                for (const auto& spawn : current_vox_spawns) {
                    Vector3 p = spawn.position;
                    // Relaxed bounds check or same strict check
                    if (p.x >= room.position.x && p.x < room.position.x + room.size.x &&
                        p.y >= room.position.y && p.y < room.position.y + room.size.y &&
                        p.z >= room.position.z && p.z < room.position.z + room.size.z) {
                        calculated_end_position = p;
                        break; 
                    }
                }
            }
            UtilityFunctions::print("LevelDensityGrid: End Position set to ", calculated_end_position, " (Room: ", room.id, ")");
        }
    }
}

void LevelDensityGrid::_apply_noise() {
     if (get_grid_size_x() <= 0 || !noise_generator.is_valid()) return;

    UtilityFunctions::print("Applying noise... Scale: ", noise_scale, ", Intensity: ", noise_intensity);
    
    // Config Updates
    noise_generator->set_seed(noise_seed);
    noise_generator->set_frequency(noise_frequency);
    
    float inv_noise_scale = (noise_scale > 0.0001) ? 1.0 / noise_scale : 1.0;
    
    int gsx = get_grid_size_x();
    int gsy = get_grid_size_y();
    int gsz = get_grid_size_z();
    float surf_thresh = get_surface_threshold();
    
    for (int z = 0; z < gsz; ++z) {
        for (int y = 0; y < gsy; ++y) {
            for (int x = 0; x < gsx; ++x) {
                Vector3i pos(x, y, z);
                float noise_val = noise_generator->get_noise_3d(static_cast<double>(x) * inv_noise_scale, static_cast<double>(y) * inv_noise_scale, static_cast<double>(z) * inv_noise_scale);
                float current_cell_value = get_cell(pos, WORLD_SOLID_VALUE);
                float new_density;
                if (current_cell_value < surf_thresh) {
                    float carving_noise = Math::min(noise_val, 0.0f);
                    new_density = current_cell_value + carving_noise * noise_intensity;
                } else {
                    new_density = current_cell_value + noise_val * noise_intensity;
                }
                set_cell(pos, new_density);
            }
        }
    }
}

// --- Vox Integration Methods ---

const ogt_vox_scene* LevelDensityGrid::_load_vox(const String &path) {
    if (path.is_empty()) return nullptr;
    
    // Check Cache
    if (vox_cache.find(path) != vox_cache.end()) {
        return vox_cache[path];
    }
    
    // Load File
    // Note: Godot's FileAccess is resource-aware (res://)
    Ref<FileAccess> file = FileAccess::open(path, FileAccess::READ);
    if (file.is_null()) {
        UtilityFunctions::printerr("LevelDensityGrid: Failed to open vox file: ", path);
        return nullptr;
    }
    
    uint64_t len = file->get_length();
    PackedByteArray buffer = file->get_buffer(len);
    
    // Parse using ogt_vox
    // We need to pass a pointer to the data. 
    // buffer.ptr() returns a const uint8_t* which is perfect.
    const ogt_vox_scene* scene = ogt_vox_read_scene(buffer.ptr(), (uint32_t)len);
    
    if (!scene) {
        UtilityFunctions::printerr("LevelDensityGrid: Failed to parse vox file: ", path);
        return nullptr;
    }
    
    // Cache and Return
    vox_cache[path] = scene;
    UtilityFunctions::print("LevelDensityGrid: Loaded vox file: ", path, " (Instances: ", scene->num_instances, ")");
    return scene;
}

void LevelDensityGrid::_carve_box(const Vector3i &origin, const Vector3i &size, float value) {
    int gsx = get_grid_size_x();
    int gsy = get_grid_size_y();
    int gsz = get_grid_size_z();
    
    // Calculate intersection
    int start_x = Math::max(0, origin.x);
    int start_y = Math::max(0, origin.y);
    int start_z = Math::max(0, origin.z);

    int end_x = Math::min(gsx, origin.x + size.x);
    int end_y = Math::min(gsy, origin.y + size.y);
    int end_z = Math::min(gsz, origin.z + size.z);

    // Safety
    if (start_x >= end_x || start_y >= end_y || start_z >= end_z) return;

    for (int z = start_z; z < end_z; ++z) {
        for (int y = start_y; y < end_y; ++y) {
            for (int x = start_x; x < end_x; ++x) {
                set_cell(Vector3i(x, y, z), value);
            }
        }
    }
}

void LevelDensityGrid::_get_vox_bounds(const ogt_vox_scene* scene, Vector3i &r_min, Vector3i &r_max) {
    r_min = Vector3i(0,0,0);
    r_max = Vector3i(0,0,0);
    if (!scene) return;
    
    int min_x = 999999, min_y = 999999, min_z = 999999;
    int max_x = -999999, max_y = -999999, max_z = -999999;
    bool has_content = false;

    for (uint32_t i = 0; i < scene->num_instances; ++i) {
        const ogt_vox_instance& inst = scene->instances[i];
        if (inst.hidden) continue;
        const ogt_vox_model* model = scene->models[inst.model_index];
        if (!model) continue;
        
        has_content = true;
        
        // Swizzle: MV X -> Godot X, MV Z -> Godot Y, MV Y -> Godot Z
        int off_x = (int)inst.transform.m30;
        int off_y_mv = (int)inst.transform.m31; 
        int off_z_mv = (int)inst.transform.m32; 
        
        int size_x = model->size_x;
        int size_y_mv = model->size_y; 
        int size_z_mv = model->size_z; 
        
        // Godot Bounds
        int start_x = off_x;           
        int end_x = off_x + size_x;
        
        int start_y = off_z_mv;          
        int end_y = off_z_mv + size_z_mv;
        
        int start_z = off_y_mv;
        int end_z = off_y_mv + size_y_mv;
        
        if (start_x < min_x) min_x = start_x;
        if (end_x > max_x) max_x = end_x;
        
        if (start_y < min_y) min_y = start_y;
        if (end_y > max_y) max_y = end_y;
        
        if (start_z < min_z) min_z = start_z;
        if (end_z > max_z) max_z = end_z;
    }
    
    if (has_content) {
        r_min = Vector3i(min_x, min_y, min_z);
        r_max = Vector3i(max_x, max_y, max_z);
    }
}

void LevelDensityGrid::_stamp_vox(const ResolvedRoom &room, const ogt_vox_scene* scene) {
    if (!scene) return;
    
    Vector3i room_origin = room.position;
    int gsx = get_grid_size_x();
    int gsy = get_grid_size_y();
    int gsz = get_grid_size_z();
    
    // Mode Logic
    if (!vox_inverse_density) {
        Vector3i min_b, max_b;
        _get_vox_bounds(scene, min_b, max_b);
        Vector3i size = max_b - min_b; // Size is Max - Min
        
        // Carve box expects Origin and Size.
        // Origin logic:
        // room_origin is the base position.
        // The Vox content exists from [min_b, max_b] relative to room_origin.
        // So we should carve from room_origin + min_b.
        
        Vector3i carve_origin = room_origin + min_b;
        
        UtilityFunctions::print("DEBUG: [Additive] Pre-carving Vox Room at ", carve_origin, " (RoomOrg ", room_origin, " + Min ", min_b, ") with Size ", size);
        _carve_box(carve_origin, size, WORLD_OPEN_VALUE);
    } else {
         UtilityFunctions::print("DEBUG: [Subtractive] Skipping Pre-carve for Vox Room at ", room_origin);
    }
    
    
    UtilityFunctions::print("DEBUG: _stamp_vox started for Room at ", room_origin, ". Instances: ", scene->num_instances);

    // Iterate Instances
    for (uint32_t i = 0; i < scene->num_instances; ++i) {
        const ogt_vox_instance& inst = scene->instances[i];
        
        // Skip hidden instances?
        if (inst.hidden) continue;

        // Get Model
        const ogt_vox_model* model = scene->models[inst.model_index];
        if (!model) continue;
        
        // Transform Logic:
        // ogt_vox flattens transforms relative to World (if read with default flags).
        // The instance transform is in 'transform.m30, m31, m32' (Translation).
        // But we need to map this to our Grid Space + Room Origin.
        
        int offset_x = (int)inst.transform.m30;
        int offset_y = (int)inst.transform.m31; // Z in Vox is Up, Y is Depth? No, Vox is Z-up usually.
        int offset_z = (int)inst.transform.m32; 
        
        // MagicaVoxel coordinate system:
        // Z is UP in MagicaVoxel. Y is Depth. X is Width.
        // Godot: Y is UP. Z is Depth. X is Width.
        // We might need to swap Y and Z if we want literal translation, 
        // OR rely on the artist exporting with "Z-Up" in mind.
        // Let's assume direct mapping for now: Vox X -> Godot X, Vox Z -> Godot Y, Vox Y -> Godot Z.
        // Wait, standard is usually:
        // MV Z -> Godot Y
        // MV Y -> Godot -Z (or Z)
        
        // Let's stick to a simple mapping first: X->X, Y->Z, Z->Y (Swizzle)
        // BUT `ogt_vox` reads raw data. 
        // Let's do a direct 1:1 mapping (X->X, Y->Y, Z->Z) first to see what happens.
        // If the user painted "Floor" on XY plane in MV, it will appear as a vertical wall in Godot (since Godot Ground is XZ).
        
        // Let's Swizzle:
        // Calculate instance bounds in grid space
        int start_grid_x = room_origin.x + offset_x;
        int start_grid_y = room_origin.y + offset_z; // Swizzle Z->Y
        int start_grid_z = room_origin.z + offset_y; // Swizzle Y->Z

        // Determine loop limits relative to local model space (0 to model->size)
        // We want: 0 <= vx < model->size_x
        // AND: 0 <= start_grid_x + vx < gsx
        // So:  -start_grid_x <= vx < gsx - start_grid_x
        
        int limit_min_vx = Math::max(0, -start_grid_x);
        int limit_max_vx = Math::min((int)model->size_x, gsx - start_grid_x);
        
        int limit_min_vy = Math::max(0, -start_grid_z); // Swizzle: vy corresponds to Z in Godot
        int limit_max_vy = Math::min((int)model->size_y, gsz - start_grid_z); // Swizzle: Model Y -> Godot Z

        int limit_min_vz = Math::max(0, -start_grid_y); // Swizzle: vz corresponds to Y in Godot
        int limit_max_vz = Math::min((int)model->size_z, gsy - start_grid_y); // Swizzle: Model Z -> Godot Y
        
        if (limit_min_vx >= limit_max_vx || limit_min_vy >= limit_max_vy || limit_min_vz >= limit_max_vz) continue;

        UtilityFunctions::print("DEBUG: Processing Instance ", i, " with clamped bounds for Grid.");
 
        for (int vx = limit_min_vx; vx < limit_max_vx; ++vx) {
            for (int vy = limit_min_vy; vy < limit_max_vy; ++vy) {
                for (int vz = limit_min_vz; vz < limit_max_vz; ++vz) {
                    
                    uint32_t voxel_index = vx + (vy * model->size_x) + (vz * model->size_x * model->size_y);
                    uint8_t color_index = model->voxel_data[voxel_index];
                    
                    // Apply swizzled coordinates
                    int final_x = start_grid_x + vx;
                    int final_y = start_grid_y + vz; // Model Z -> Godot Y
                    int final_z = start_grid_z + vy; // Model Y -> Godot Z
                    
                    Vector3i grid_pos(final_x, final_y, final_z);

                    if (color_index == 0) { // MV AIR
                        if (!vox_inverse_density) {
                            // Additive: Explicitly ensure it is Air
                            set_cell(grid_pos, WORLD_OPEN_VALUE);
                        }
                    } else { // MV SOLID
                        if (current_vox_spawn_map.has(color_index)) { // Spawn Point
                            set_cell(grid_pos, WORLD_OPEN_VALUE);
                            VoxSpawn spawn;
                            spawn.position = Vector3(final_x, final_y, final_z);
                            spawn.palette_index = color_index;
                            current_vox_spawns.push_back(spawn);
                        } else { // Block
                            if (vox_inverse_density) {
                                // Subtractive
                                set_cell(grid_pos, WORLD_OPEN_VALUE);
                            } else {
                                // Additive
                                set_cell(grid_pos, WORLD_SOLID_VALUE);
                                if (current_vox_material_map.has(color_index)) {
                                    int mat_id = current_vox_material_map[color_index];
                                    set_material_id(grid_pos, mat_id);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void LevelDensityGrid::_clear_vox_cache() {
    for (auto const& [path, scene] : vox_cache) {
        ogt_vox_destroy_scene(scene);
    }
    vox_cache.clear();
}

TypedArray<Dictionary> LevelDensityGrid::get_vox_spawns() const {
    TypedArray<Dictionary> ret;
    for (const auto& s : current_vox_spawns) {
        Dictionary d;
        d["position"] = s.position;
        d["palette_index"] = s.palette_index;
        ret.push_back(d);
    }
    return ret;
}



void LevelDensityGrid::low_pass() {
    int smoothing = smoothing_strength;
    int gsx = get_grid_size_x();
    int gsy = get_grid_size_y();
    int gsz = get_grid_size_z();
    int required_size = smoothing * 2 + 1;
    if (gsx < required_size || gsy < required_size || gsz < required_size) {
        UtilityFunctions::printerr("LevelDensityGrid.low_pass: Grid is too small for the given smoothing strength.");
        return;
    }
    std::vector<float> smoothed_grid;
    size_t grid_total_size = (size_t)gsx * gsy * gsz;
    
    smoothed_grid.resize(grid_total_size);
    for (int z = smoothing; z < gsz - smoothing; ++z) {
        for (int y = smoothing; y < gsy - smoothing; ++y) {
            for (int x = smoothing; x < gsx - smoothing; ++x) {
                Vector3i current_pos(x, y, z);
                float sum = 0.0f;
                int count = 0;
                for (int dz = -smoothing; dz <= smoothing; ++dz) {
                    for (int dy = -smoothing; dy <= smoothing; ++dy) {
                        for (int dx = -smoothing; dx <= smoothing; ++dx) {
                            Vector3i neighbor_pos = current_pos + Vector3i(dx, dy, dz);
                            sum += get_cell(neighbor_pos, WORLD_SOLID_VALUE);
                            count++;
                        }
                    }
                }
                float average = (count > 0) ? sum / count : get_cell(current_pos, WORLD_SOLID_VALUE);
                int index = get_index(current_pos);
                 if (index != -1) {
                    smoothed_grid[index] = average;
                 }
            }
        }
    }
    for (int z = smoothing; z < gsz - smoothing; ++z) {
        for (int y = smoothing; y < gsy - smoothing; ++y) {
            for (int x = smoothing; x < gsx - smoothing; ++x) {
                Vector3i pos(x, y, z);
                int index = get_index(pos);
                 if (index != -1) {
                     set_cell(pos, smoothed_grid[index]);
                 }
            }
        }
    }
}

void LevelDensityGrid::_calculate_surface_normals() {
    surface_normals.clear();
    int gsx = get_grid_size_x();
    int gsy = get_grid_size_y();
    int gsz = get_grid_size_z();
    float surf_thresh = get_surface_threshold();

    UtilityFunctions::print("Calculating surface normals...");
    
    // (Existing logic preserved)
    for (int x = 0; x < gsx; ++x) {
        for (int y = 0; y < gsy; ++y) {
            for (int z = 0; z < gsz; ++z) {
                Vector3i current_pos(x, y, z);
                if (get_cell(current_pos, WORLD_SOLID_VALUE) > surf_thresh) {
                    bool is_surface = false;
                    Vector3 gradient;
                    if (x > 0) {
                        float left = get_cell(Vector3i(x - 1, y, z), WORLD_SOLID_VALUE);
                        if (left < surf_thresh) is_surface = true;
                        gradient.x -= left;
                    }
                    if (x < gsx - 1) {
                        float right = get_cell(Vector3i(x + 1, y, z), WORLD_SOLID_VALUE);
                        if (right < surf_thresh) is_surface = true;
                        gradient.x += right;
                    }
                    // Check Y neighbors...
                    if (y > 0) {
                        float down = get_cell(Vector3i(x, y - 1, z), WORLD_SOLID_VALUE);
                        if (down < surf_thresh) is_surface = true;
                        gradient.y -= down;
                    }
                    if (y < gsy - 1) {
                        float up = get_cell(Vector3i(x, y + 1, z), WORLD_SOLID_VALUE);
                        if (up < surf_thresh) is_surface = true;
                        gradient.y += up;
                    }
                    // Check Z neighbors...
                    if (z > 0) {
                        float back = get_cell(Vector3i(x, y, z - 1), WORLD_SOLID_VALUE);
                        if (back < surf_thresh) is_surface = true;
                        gradient.z -= back;
                    }
                    if (z < gsz - 1) {
                        float front = get_cell(Vector3i(x, y, z + 1), WORLD_SOLID_VALUE);
                        if (front < surf_thresh) is_surface = true;
                        gradient.z += front;
                    }

                    if (is_surface) {
                        Vector3 normal = gradient.normalized();
                        surface_normals[current_pos] = -normal; // Store normal
                    }
                }
            }
        }
    }
}

void LevelDensityGrid::_mark_brush(const Vector3i &center, int radius_lower, int radius_higher, float value) {
    // Forward to PathCarver's logic? Or duplicate?
    // PathCarver has a private `_mark_brush`.
    // Since _mark_brush is exposed to GDScript, we should implement it here for manual usage.
    // We can copy the logic or friend `PathCarver`.
    // Copying simple logic is fine to avoid tight coupling hacks.
    
    // We reuse the logic from PathCarver but adapted for this class scope
    // Actually, can we just use `path_carver` instance?
    // `path_carver._mark_brush` is private.
    // Let's reimplement simple logic here. It's just a paint loop.
    
    int r = radius_lower; // simplified
    int gsx = get_grid_size_x();
    int gsy = get_grid_size_y();
    int gsz = get_grid_size_z();
    
    // Simple Box or Sphere
    for(int x = -r; x <= r; x++) {
        for(int y = -r; y <= r; y++) {
            for(int z = -r; z <= r; z++) {
                Vector3i p = center + Vector3i(x,y,z);
                if (path_carver.use_square_brush) { // reuse property
                     if (p.x >= 0 && p.x < gsx && p.y >= 0 && p.y < gsy && p.z >= 0 && p.z < gsz) {
                         set_cell(p, value);
                     }
                } else {
                    if (x*x + y*y + z*z <= r*r) {
                        if (p.x >= 0 && p.x < gsx && p.y >= 0 && p.y < gsy && p.z >= 0 && p.z < gsz) {
                            set_cell(p, value);
                        }
                    }
                }
            }
        }
    }
}

Vector3i LevelDensityGrid::_find_ground_position(const Vector3i &start_pos) {
    // Forward to PathCarver logic? Or just same loop.
     int gsx = get_grid_size_x();
    int gsy = get_grid_size_y();
    int gsz = get_grid_size_z();
    
    int x = Math::clamp(start_pos.x, 0, gsx - 1);
    int z = Math::clamp(start_pos.z, 0, gsz - 1);
    int start_y = Math::clamp(start_pos.y, 0, gsy - 1);

    for (int y = start_y; y >= 0; --y) {
        Vector3i pos(x, y, z);
        if (get_cell(pos, WORLD_SOLID_VALUE) >= get_surface_threshold()) {
            return Vector3i(x, y + 1, z);
        }
    }
    return Vector3i(x, 0, z);
}

Ref<DensityGrid> LevelDensityGrid::generate_liquid_grid() {
    return liquid_gen.generate_liquid_grid(this);
}



// Noise
void LevelDensityGrid::set_noise_scale(float p_scale) { noise_scale = p_scale; }
float LevelDensityGrid::get_noise_scale() const { return noise_scale; }
void LevelDensityGrid::set_noise_intensity(float p_intensity) { noise_intensity = p_intensity; }
float LevelDensityGrid::get_noise_intensity() const { return noise_intensity; }
void LevelDensityGrid::set_noise_seed(int p_seed) {
    noise_seed = p_seed;
    if (noise_generator.is_valid()) noise_generator->set_seed(noise_seed);
    // Also update liquid gen dependency?
    liquid_gen.noise_seed = p_seed; 
}
int LevelDensityGrid::get_noise_seed() const { return noise_seed; }
void LevelDensityGrid::set_noise_frequency(float p_freq) {
    noise_frequency = p_freq;
     if (noise_generator.is_valid()) noise_generator->set_frequency(noise_frequency);
}
float LevelDensityGrid::get_noise_frequency() const { return noise_frequency; }
void LevelDensityGrid::set_noise_generator(const Ref<FastNoiseLite> &p_noise) {
    noise_generator = p_noise;
    if (noise_generator.is_valid()) {
        noise_generator->set_seed(noise_seed);
        noise_generator->set_frequency(noise_frequency);
    }
}
Ref<FastNoiseLite> LevelDensityGrid::get_noise_generator() const { return noise_generator; }

// Room Gen
void LevelDensityGrid::set_room_count(int p_count) { room_gen.room_count = p_count > 0 ? p_count : 0; }
int LevelDensityGrid::get_room_count() const { return room_gen.room_count; }
void LevelDensityGrid::set_min_room_size(const Vector3i &p_size) { room_gen.min_room_size = p_size; }
Vector3i LevelDensityGrid::get_min_room_size() const { return room_gen.min_room_size; }
void LevelDensityGrid::set_max_room_size(const Vector3i &p_size) { room_gen.max_room_size = p_size; }
Vector3i LevelDensityGrid::get_max_room_size() const { return room_gen.max_room_size; }
void LevelDensityGrid::set_max_placement_tries(int p_tries) { room_gen.max_placement_tries = p_tries > 0 ? p_tries : 1; }
int LevelDensityGrid::get_max_placement_tries() const { return room_gen.max_placement_tries; }

// Path Carver
void LevelDensityGrid::set_path_brush_min_radius(int p_radius) { path_carver.path_brush_min_radius = p_radius >= 0 ? p_radius : 0; }
int LevelDensityGrid::get_path_brush_min_radius() const { return path_carver.path_brush_min_radius; }
void LevelDensityGrid::set_path_brush_max_radius(int p_radius) { path_carver.path_brush_max_radius = p_radius >= 0 ? p_radius : 0; }
int LevelDensityGrid::get_path_brush_max_radius() const { return path_carver.path_brush_max_radius; }
void LevelDensityGrid::set_use_square_brush(bool p_enabled) { path_carver.use_square_brush = p_enabled; }
bool LevelDensityGrid::get_use_square_brush() const { return path_carver.use_square_brush; }
void LevelDensityGrid::set_vertical_movement_cost_multiplier(float p_mult) { path_carver.vertical_movement_cost_multiplier = p_mult > 0 ? p_mult : 1.0; }
float LevelDensityGrid::get_vertical_movement_cost_multiplier() const { return path_carver.vertical_movement_cost_multiplier; }
void LevelDensityGrid::set_dungeon_mode(bool p_enabled) { path_carver.dungeon_mode = p_enabled; }
bool LevelDensityGrid::get_dungeon_mode() const { return path_carver.dungeon_mode; }
void LevelDensityGrid::set_dungeon_path_algorithm(int p_algo) { path_carver.dungeon_path_algorithm = p_algo; }
int LevelDensityGrid::get_dungeon_path_algorithm() const { return path_carver.dungeon_path_algorithm; }

// Bezier
void LevelDensityGrid::set_path_segments(int p_segments) { path_carver.path_segments = p_segments > 0 ? p_segments : 1; }
int LevelDensityGrid::get_path_segments() const { return path_carver.path_segments; }

void LevelDensityGrid::set_path_bend_factor(float p_factor) { path_carver.path_bend_factor = p_factor >= 0.0f ? p_factor : 0.0f; }
float LevelDensityGrid::get_path_bend_factor() const { return path_carver.path_bend_factor; }

void LevelDensityGrid::set_path_wobble_magnitude(float p_mag) { path_carver.path_wobble_magnitude = p_mag >= 0.0f ? p_mag : 0.0f; }
float LevelDensityGrid::get_path_wobble_magnitude() const { return path_carver.path_wobble_magnitude; }

void LevelDensityGrid::set_path_wobble_frequency(float p_freq) { path_carver.path_wobble_frequency = p_freq; }
float LevelDensityGrid::get_path_wobble_frequency() const { return path_carver.path_wobble_frequency; }

void LevelDensityGrid::set_connect_from_ground_level(bool p_enabled) { path_carver.connect_from_ground_level = p_enabled; }
bool LevelDensityGrid::get_connect_from_ground_level() const { return path_carver.connect_from_ground_level; }



// Smoothing (On Helper)
void LevelDensityGrid::set_smooth_terrain(bool p_enabled) { smooth_terrain = p_enabled; }
bool LevelDensityGrid::get_smooth_terrain() const { return smooth_terrain; }
void LevelDensityGrid::set_smoothing_strength(int p_strength) { smoothing_strength = p_strength > 0 ? p_strength : 1; }
int LevelDensityGrid::get_smoothing_strength() const { return smoothing_strength; }

// Liquid Gen
void LevelDensityGrid::set_max_basin_size(int p_size) { liquid_gen.max_basin_size = p_size > 1 ? p_size : 1; }
int LevelDensityGrid::get_max_basin_size() const { return liquid_gen.max_basin_size; }
void LevelDensityGrid::set_sparsity_cutoff(float p_cutoff) { liquid_gen.sparsity_cutoff = p_cutoff; }
float LevelDensityGrid::get_sparsity_cutoff() const { return liquid_gen.sparsity_cutoff; }
void LevelDensityGrid::set_water_height_density(float p_density) { liquid_gen.water_height_density = Math::clamp(p_density, 0.0f, 1.0f); }
float LevelDensityGrid::get_water_height_density() const { return liquid_gen.water_height_density; }
void LevelDensityGrid::set_liquid_resolution_multiplier(int p_mult) { liquid_gen.liquid_resolution_multiplier = (p_mult < 1) ? 1 : p_mult; }
int LevelDensityGrid::get_liquid_resolution_multiplier() const { return liquid_gen.liquid_resolution_multiplier; }

// Runtime
Vector3 LevelDensityGrid::get_calculated_spawn_position() const { return calculated_spawn_position; }
Vector3 LevelDensityGrid::get_calculated_end_position() const { return calculated_end_position; }
Dictionary LevelDensityGrid::get_surface_normals() const { return surface_normals; }

String LevelDensityGrid::get_room_type_at(const Vector3 &world_pos, float voxel_size) const {
    Vector3i grid_pos = Vector3i(world_pos / voxel_size);
    for(const auto &room : resolved_rooms) {
        if (grid_pos.x >= room.position.x && grid_pos.x < room.position.x + room.size.x &&
            grid_pos.y >= room.position.y && grid_pos.y < room.position.y + room.size.y &&
            grid_pos.z >= room.position.z && grid_pos.z < room.position.z + room.size.z) {
            return room.type;
        }
    }
    return "corridor"; 
}



// Vox Integration
void LevelDensityGrid::set_vox_inverse_density(bool p_enabled) { vox_inverse_density = p_enabled; }
bool LevelDensityGrid::get_vox_inverse_density() const { return vox_inverse_density; }

void LevelDensityGrid::_bind_methods() {
    // --- Main Generation Function ---
    ClassDB::bind_method(D_METHOD("generate_level_data", "world_grid_dimensions", "voxel_size", "seed"), &LevelDensityGrid::generate_level_data);
    ClassDB::bind_method(D_METHOD("generate_from_graph", "node_data", "edge_data", "voxel_size", "seed", "settings"), &LevelDensityGrid::generate_from_graph, DEFVAL(Dictionary()));
    ClassDB::bind_method(D_METHOD("_mark_brush", "center", "radius_lower", "radius_higher", "value"), &LevelDensityGrid::_mark_brush);
    // --- Getters for Runtime Data ---
    ClassDB::bind_method(D_METHOD("get_calculated_spawn_position"), &LevelDensityGrid::get_calculated_spawn_position);
    ClassDB::bind_method(D_METHOD("get_calculated_end_position"), &LevelDensityGrid::get_calculated_end_position);
    ClassDB::bind_method(D_METHOD("get_surface_normals"), &LevelDensityGrid::get_surface_normals);
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3, "calculated_spawn_position", PROPERTY_HINT_NONE, "", PROPERTY_USAGE_READ_ONLY), "", "get_calculated_spawn_position");
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3, "calculated_end_position", PROPERTY_HINT_NONE, "", PROPERTY_USAGE_READ_ONLY), "", "get_calculated_end_position");
    ADD_PROPERTY(PropertyInfo(Variant::DICTIONARY, "surface_normals", PROPERTY_HINT_NONE, "", PROPERTY_USAGE_READ_ONLY), "", "get_surface_normals");


    // --- Property Getters/Setters ---
    // Noise Group
    ADD_GROUP("Noise", "noise_");
    ClassDB::bind_method(D_METHOD("set_noise_scale", "scale"), &LevelDensityGrid::set_noise_scale);
    ClassDB::bind_method(D_METHOD("get_noise_scale"), &LevelDensityGrid::get_noise_scale);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "noise_scale"), "set_noise_scale", "get_noise_scale");

    ClassDB::bind_method(D_METHOD("set_noise_intensity", "intensity"), &LevelDensityGrid::set_noise_intensity);
    ClassDB::bind_method(D_METHOD("get_noise_intensity"), &LevelDensityGrid::get_noise_intensity);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "noise_intensity"), "set_noise_intensity", "get_noise_intensity");

    ClassDB::bind_method(D_METHOD("set_noise_seed", "seed"), &LevelDensityGrid::set_noise_seed);
    ClassDB::bind_method(D_METHOD("get_noise_seed"), &LevelDensityGrid::get_noise_seed);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "noise_seed"), "set_noise_seed", "get_noise_seed");

    ClassDB::bind_method(D_METHOD("set_noise_frequency", "freq"), &LevelDensityGrid::set_noise_frequency);
    ClassDB::bind_method(D_METHOD("get_noise_frequency"), &LevelDensityGrid::get_noise_frequency);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "noise_frequency"), "set_noise_frequency", "get_noise_frequency");

    ClassDB::bind_method(D_METHOD("set_noise_generator", "noise"), &LevelDensityGrid::set_noise_generator);
    ClassDB::bind_method(D_METHOD("get_noise_generator"), &LevelDensityGrid::get_noise_generator);
    ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "noise_generator", PROPERTY_HINT_RESOURCE_TYPE, "FastNoiseLite"), "set_noise_generator", "get_noise_generator");

    // Rooms Group
    ADD_GROUP("Rooms", "room_");
    ClassDB::bind_method(D_METHOD("set_room_count", "count"), &LevelDensityGrid::set_room_count);
    ClassDB::bind_method(D_METHOD("get_room_count"), &LevelDensityGrid::get_room_count);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "room_count"), "set_room_count", "get_room_count");

    ClassDB::bind_method(D_METHOD("set_min_room_size", "size"), &LevelDensityGrid::set_min_room_size);
    ClassDB::bind_method(D_METHOD("get_min_room_size"), &LevelDensityGrid::get_min_room_size);
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3I, "room_min_size"), "set_min_room_size", "get_min_room_size");

    ClassDB::bind_method(D_METHOD("set_max_room_size", "size"), &LevelDensityGrid::set_max_room_size);
    ClassDB::bind_method(D_METHOD("get_max_room_size"), &LevelDensityGrid::get_max_room_size);
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3I, "room_max_size"), "set_max_room_size", "get_max_room_size");

    ClassDB::bind_method(D_METHOD("set_max_placement_tries", "tries"), &LevelDensityGrid::set_max_placement_tries);
    ClassDB::bind_method(D_METHOD("get_max_placement_tries"), &LevelDensityGrid::get_max_placement_tries);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "room_max_placement_tries"), "set_max_placement_tries", "get_max_placement_tries");


    // Paths Group - CLEANED UP WITH SUBGROUPS
    ADD_GROUP("Paths", "path_");
    
    // General Path Settings
    ClassDB::bind_method(D_METHOD("set_dungeon_mode", "enabled"), &LevelDensityGrid::set_dungeon_mode);
    ClassDB::bind_method(D_METHOD("get_dungeon_mode"), &LevelDensityGrid::get_dungeon_mode);
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "path_dungeon_mode"), "set_dungeon_mode", "get_dungeon_mode");

    ClassDB::bind_method(D_METHOD("set_path_brush_min_radius", "radius"), &LevelDensityGrid::set_path_brush_min_radius);
    ClassDB::bind_method(D_METHOD("get_path_brush_min_radius"), &LevelDensityGrid::get_path_brush_min_radius);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "path_brush_min_radius"), "set_path_brush_min_radius", "get_path_brush_min_radius");

    ClassDB::bind_method(D_METHOD("set_path_brush_max_radius", "radius"), &LevelDensityGrid::set_path_brush_max_radius);
    ClassDB::bind_method(D_METHOD("get_path_brush_max_radius"), &LevelDensityGrid::get_path_brush_max_radius);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "path_brush_max_radius"), "set_path_brush_max_radius", "get_path_brush_max_radius");

    ClassDB::bind_method(D_METHOD("set_use_square_brush", "enabled"), &LevelDensityGrid::set_use_square_brush);
    ClassDB::bind_method(D_METHOD("get_use_square_brush"), &LevelDensityGrid::get_use_square_brush);
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "path_use_square_brush"), "set_use_square_brush", "get_use_square_brush");

    ClassDB::bind_method(D_METHOD("set_vertical_movement_cost_multiplier", "mult"), &LevelDensityGrid::set_vertical_movement_cost_multiplier);
    ClassDB::bind_method(D_METHOD("get_vertical_movement_cost_multiplier"), &LevelDensityGrid::get_vertical_movement_cost_multiplier);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "path_vertical_movement_cost_multiplier"), "set_vertical_movement_cost_multiplier", "get_vertical_movement_cost_multiplier");
    



    // Dungeon Specific
    ADD_SUBGROUP("Dungeon Settings", "path_dungeon_");
    ClassDB::bind_method(D_METHOD("set_dungeon_path_algorithm", "algo"), &LevelDensityGrid::set_dungeon_path_algorithm);
    ClassDB::bind_method(D_METHOD("get_dungeon_path_algorithm"), &LevelDensityGrid::get_dungeon_path_algorithm);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "dungeon_path_algorithm", PROPERTY_HINT_ENUM, "AStar (Slow/Optimal),Castle (Fast/Recursive)"), "set_dungeon_path_algorithm", "get_dungeon_path_algorithm");

    // Bezier Specific
    ADD_SUBGROUP("Bezier Settings", "path_bezier_");
    ClassDB::bind_method(D_METHOD("set_path_segments", "segments"), &LevelDensityGrid::set_path_segments);
    ClassDB::bind_method(D_METHOD("get_path_segments"), &LevelDensityGrid::get_path_segments);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "path_segments", PROPERTY_HINT_RANGE, "1,10,1"), "set_path_segments", "get_path_segments");

    ClassDB::bind_method(D_METHOD("set_path_bend_factor", "factor"), &LevelDensityGrid::set_path_bend_factor);
    ClassDB::bind_method(D_METHOD("get_path_bend_factor"), &LevelDensityGrid::get_path_bend_factor);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "path_bend_factor", PROPERTY_HINT_RANGE, "0.0,2.0,0.05"), "set_path_bend_factor", "get_path_bend_factor");
    
    ClassDB::bind_method(D_METHOD("set_path_wobble_magnitude", "magnitude"), &LevelDensityGrid::set_path_wobble_magnitude);
    ClassDB::bind_method(D_METHOD("get_path_wobble_magnitude"), &LevelDensityGrid::get_path_wobble_magnitude);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "path_wobble_magnitude", PROPERTY_HINT_RANGE, "0.0,10.0,0.1"), "set_path_wobble_magnitude", "get_path_wobble_magnitude");
    
    ClassDB::bind_method(D_METHOD("set_path_wobble_frequency", "frequency"), &LevelDensityGrid::set_path_wobble_frequency);
    ClassDB::bind_method(D_METHOD("get_path_wobble_frequency"), &LevelDensityGrid::get_path_wobble_frequency);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "path_wobble_frequency", PROPERTY_HINT_RANGE, "0.0,1.0,0.01"), "set_path_wobble_frequency", "get_path_wobble_frequency");

    ClassDB::bind_method(D_METHOD("set_connect_from_ground_level", "enabled"), &LevelDensityGrid::set_connect_from_ground_level);
    ClassDB::bind_method(D_METHOD("get_connect_from_ground_level"), &LevelDensityGrid::get_connect_from_ground_level);
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "path_connect_from_ground_level"), "set_connect_from_ground_level", "get_connect_from_ground_level");


    // Smoothing Group
    ADD_GROUP("Smoothing", "smoothing_");
    ClassDB::bind_method(D_METHOD("set_smooth_terrain", "enabled"), &LevelDensityGrid::set_smooth_terrain);
    ClassDB::bind_method(D_METHOD("get_smooth_terrain"), &LevelDensityGrid::get_smooth_terrain);
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "smoothing_enabled"), "set_smooth_terrain", "get_smooth_terrain");
    
    ClassDB::bind_method(D_METHOD("set_smoothing_strength", "strength"), &LevelDensityGrid::set_smoothing_strength);
    ClassDB::bind_method(D_METHOD("get_smoothing_strength"), &LevelDensityGrid::get_smoothing_strength);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "smoothing_strength", PROPERTY_HINT_RANGE, "1,5,1"), "set_smoothing_strength", "get_smoothing_strength");


    ClassDB::bind_method(D_METHOD("_find_ground_position", "start_pos"), &LevelDensityGrid::_find_ground_position);
    ClassDB::bind_method(D_METHOD("generate_liquid_grid"), &LevelDensityGrid::generate_liquid_grid);
    
    ADD_GROUP("Liquid Generation", "liquid_");
    
    ClassDB::bind_method(D_METHOD("set_max_basin_size", "size"), &LevelDensityGrid::set_max_basin_size);
    ClassDB::bind_method(D_METHOD("get_max_basin_size"), &LevelDensityGrid::get_max_basin_size);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "liquid_max_basin_size", PROPERTY_HINT_RANGE, "1,500,1"), "set_max_basin_size", "get_max_basin_size");

    ClassDB::bind_method(D_METHOD("set_sparsity_cutoff", "cutoff"), &LevelDensityGrid::set_sparsity_cutoff);
    ClassDB::bind_method(D_METHOD("get_sparsity_cutoff"), &LevelDensityGrid::get_sparsity_cutoff);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "liquid_sparsity_cutoff", PROPERTY_HINT_RANGE, "-1.0,1.0,0.05"), "set_sparsity_cutoff", "get_sparsity_cutoff");

    // VOX
    ClassDB::bind_method(D_METHOD("get_vox_spawns"), &LevelDensityGrid::get_vox_spawns);
    ClassDB::bind_method(D_METHOD("set_vox_inverse_density", "enabled"), &LevelDensityGrid::set_vox_inverse_density);
    ClassDB::bind_method(D_METHOD("get_vox_inverse_density"), &LevelDensityGrid::get_vox_inverse_density);
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "vox_inverse_density"), "set_vox_inverse_density", "get_vox_inverse_density");

    ClassDB::bind_method(D_METHOD("set_placement_algorithm", "algorithm"), &LevelDensityGrid::set_placement_algorithm);
    ClassDB::bind_method(D_METHOD("get_placement_algorithm"), &LevelDensityGrid::get_placement_algorithm);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "placement_algorithm", PROPERTY_HINT_ENUM, "Force Directed,BSP"), "set_placement_algorithm", "get_placement_algorithm");

    BIND_ENUM_CONSTANT(ALGO_FORCE_DIRECTED);
    BIND_ENUM_CONSTANT(ALGO_BSP);

    // Debug
    ClassDB::bind_method(D_METHOD("get_surface_zones"), &LevelDensityGrid::get_surface_zones);
    ClassDB::bind_method(D_METHOD("populate_multimesh", "multimesh", "points", "scale"), &LevelDensityGrid::populate_multimesh);

    ClassDB::bind_method(D_METHOD("set_water_height_density", "density"), &LevelDensityGrid::set_water_height_density);
    ClassDB::bind_method(D_METHOD("get_water_height_density"), &LevelDensityGrid::get_water_height_density);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "liquid_water_height_density", PROPERTY_HINT_RANGE, "0.0,1.0,0.01"), "set_water_height_density", "get_water_height_density");

// ... existing code ...
    ClassDB::bind_method(D_METHOD("set_liquid_resolution_multiplier", "multiplier"), &LevelDensityGrid::set_liquid_resolution_multiplier);
    ClassDB::bind_method(D_METHOD("get_liquid_resolution_multiplier"), &LevelDensityGrid::get_liquid_resolution_multiplier);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "liquid_resolution_multiplier", PROPERTY_HINT_RANGE, "1,8,1"), "set_liquid_resolution_multiplier", "get_liquid_resolution_multiplier");


}

void LevelDensityGrid::set_placement_algorithm(PlacementAlgorithm p_algo) {
    placement_algorithm = p_algo;
}

LevelDensityGrid::PlacementAlgorithm LevelDensityGrid::get_placement_algorithm() const {
    return placement_algorithm;
}


void LevelDensityGrid::_place_rooms_force_directed(std::vector<ResolvedRoom> &processing_rooms, const std::vector<ResolvedEdge> &edges, const std::vector<int> &node_depths, int max_depth, const TypedArray<Dictionary> &macro_volumes, const PackedVector3Array &guide_points) {
    // Reconstruct id_to_index for constraints 
    // (Ideally pass this map, but for now we rebuild it quickly)
    std::map<String, int> id_to_index;
    for(int i=0; i<processing_rooms.size(); ++i) {
        id_to_index[processing_rooms[i].id] = i;
    }

    // Config
    int iterations = 100;
    float repulsion = 15.0f;
    float spring_length = 30.0f;
    float spring_strength = 0.5f;
    Vector3 center_target = Vector3(get_grid_size_x()/2, get_grid_size_y()/2, get_grid_size_z()/2);
    
    // Guide Curve Lambda
    auto sample_curve = [&](float t) -> Vector3 {
        if (guide_points.size() < 2) return Vector3(0,0,0);
        float total_len = guide_points.size() - 1.0f;
        float segment = t * total_len;
        int idx = (int)segment;
        if (idx >= guide_points.size() - 1) return guide_points[guide_points.size() - 1];
        float sub_t = segment - idx;
        return guide_points[idx].lerp(guide_points[idx+1], sub_t);
    };

    for(int k=0; k<iterations; ++k) {
        std::vector<Vector3> forces(processing_rooms.size(), Vector3(0,0,0));

        for(int i=0; i<processing_rooms.size(); ++i) {
            if (processing_rooms[i].is_fixed) continue; 

            for(int j=i+1; j<processing_rooms.size(); ++j) {
                Vector3 diff = processing_rooms[i].center() - processing_rooms[j].center();
                float dist = diff.length();
                if(dist < 0.1f) dist = 0.1f;
                float combined_size = (processing_rooms[i].size.x + processing_rooms[j].size.x) * 0.8f;
                if (dist < combined_size) {
                    float force_mag = repulsion * (combined_size / dist);
                    Vector3 f = diff.normalized() * force_mag;
                    forces[i] += f;
                    forces[j] -= f;
                }
            }
            Vector3 to_center = center_target - processing_rooms[i].center();
            to_center.y *= 0.1f; // Weak vertical gravity
            forces[i] += to_center * 0.05f; 
            
            // Relative Constraints
            if (processing_rooms[i].has_relative_constraint) {
                String target_id = processing_rooms[i].relative_target_id;
                if (id_to_index.count(target_id)) {
                    int target_idx = id_to_index[target_id];
                    Vector3 target_center = processing_rooms[target_idx].center();
                    Vector3 desired_pos = target_center + processing_rooms[i].relative_offset;
                    Vector3 current_center = processing_rooms[i].center();
                    Vector3 diff = desired_pos - current_center;
                    forces[i] += diff * 2.0f; 
                }
            }

            // Macro Volumes
            for(int v=0; v < macro_volumes.size(); ++v) {
                Dictionary vol = macro_volumes[v];
                String tag = vol.get("tag", "");
                if (tag == processing_rooms[i].type) {
                     AABB bounds = vol.get("bounds", AABB());
                     Vector3 center = processing_rooms[i].center();
                     if (!bounds.has_point(center)) {
                         Vector3 clamped;
                         clamped.x = Math::clamp(center.x, bounds.position.x, bounds.position.x + bounds.size.x);
                         clamped.y = Math::clamp(center.y, bounds.position.y, bounds.position.y + bounds.size.y);
                         clamped.z = Math::clamp(center.z, bounds.position.z, bounds.position.z + bounds.size.z);
                         forces[i] += (clamped - center) * 1.0f; 
                     }
                }
            }
        }

        for(const auto &edge : edges) {
            int idxA = edge.from_index; 
            int idxB = edge.to_index;   
            Vector3 diff = processing_rooms[idxB].center() - processing_rooms[idxA].center();
            float dist = diff.length();
            float displacement = dist - spring_length;
            Vector3 f = diff.normalized() * (displacement * spring_strength);
            forces[idxA] += f;
            forces[idxB] -= f;
        }

        for(int i=0; i<processing_rooms.size(); ++i) {
            if (processing_rooms[i].is_fixed) {
                 Vector3 fix_center = processing_rooms[i].fixed_position;
                  processing_rooms[i].position = Vector3i(
                    fix_center.x - processing_rooms[i].size.x / 2,
                    fix_center.y - processing_rooms[i].size.y / 2,
                    fix_center.z - processing_rooms[i].size.z / 2
                );
                continue;
            }

            if (forces[i].length() > 10.0f) forces[i] = forces[i].normalized() * 10.0f;
            Vector3 new_center = processing_rooms[i].center() + forces[i];
            
            // Allow bounds to be a bit looser during sim, or clamp strictly? 
            // Original code clamped strictly.
            processing_rooms[i].position = Vector3i(
                Math::clamp((int)new_center.x - processing_rooms[i].size.x / 2, 0, get_grid_size_x() - processing_rooms[i].size.x),
                Math::clamp((int)new_center.y, 2, get_grid_size_y() - 5) - processing_rooms[i].size.y / 2,
                Math::clamp((int)new_center.z - processing_rooms[i].size.z / 2, 0, get_grid_size_z() - processing_rooms[i].size.z)
            );
        }
    }
}

void LevelDensityGrid::_place_rooms_bsp(std::vector<ResolvedRoom> &processing_rooms) {
    if (processing_rooms.empty()) return;

    // 1. Initial Bounds with 1-Voxel Margin
    int gsx = get_grid_size_x();
    int gsy = get_grid_size_y();
    int gsz = get_grid_size_z();
    
    // AABB in integer grid space
    struct GridAABB {
        Vector3i min;
        Vector3i max; // Inclusive? Or Size? Let's use Min and Size.
        Vector3i size;
        
        int volume() const { return size.x * size.y * size.z; }
    };
    
    // Start with Margin: [1, 1, 1] to [gsx-1, gsy-1, gsz-1] (Exclusive max)
    // Size = BoundMax - BoundMin
    // Size = (gsx-1) - 1 = gsx - 2
    GridAABB root;
    root.min = Vector3i(1, 1, 1);
    root.size = Vector3i(gsx - 2, gsy - 2, gsz - 2);
    
    if (root.size.x <= 0 || root.size.y <= 0 || root.size.z <= 0) {
        UtilityFunctions::printerr("Grid too small for BSP with margins.");
        return;
    }

    std::vector<GridAABB> leaves;
    leaves.push_back(root);

    // 2. Recursive Split
    int rooms_needed = processing_rooms.size();
    int safety_max = 10000;
    int iterations = 0;

    while (leaves.size() < rooms_needed && iterations < safety_max) {
        iterations++;
        
        // Find largest leaf
        int best_idx = -1;
        int max_vol = -1;
        
        for (int i=0; i<leaves.size(); ++i) {
            int vol = leaves[i].volume();
            if (vol > max_vol) {
                max_vol = vol;
                best_idx = i;
            }
        }
        
        if (best_idx == -1) break; // Should not happen
        
        GridAABB current = leaves[best_idx];
        
        // Attempt split
        // Pick axis based on longest dimension
        int axis = 0; // 0=x, 1=y, 2=z
        if (current.size.y > current.size.x && current.size.y > current.size.z) axis = 1;
        else if (current.size.z > current.size.x) axis = 2;
        
        // Randomize axis occasionally for variety?
        if (rng->randf() > 0.8) axis = rng->randi() % 3;

        // Determine split point
        // Min size for a child node? Let's say 5 voxels?
        int min_child_size = 8; 
        
        int current_len = (axis == 0) ? current.size.x : (axis == 1) ? current.size.y : current.size.z;
        
        if (current_len < min_child_size * 2) {
             // Too small to split this axis. Try another? 
             // Just implicitly marks it as unsplittable by loop selection logic? 
             // Ideally we should mark it as "do not split" or prioritize valid splits.
             // For simplicity, we just stop splitting this node if selected. 
             // But if we need more rooms, we must find *some* node to split.
             // Let's degrade `max_vol` check to ignore too-small nodes.
             // Simpler approach: If selected node is too small, skip it and try next largest?
             // Since we iterate `leaves`, we modify list.
             // Let's fail-safe: just break if we can't split the biggest.
             break; 
        }

        int split_pos = rng->randi_range(min_child_size, current_len - min_child_size);
        
        // Create Children
        GridAABB c1 = current;
        GridAABB c2 = current;
        
        if (axis == 0) {
            c1.size.x = split_pos;
            c2.min.x += split_pos;
            c2.size.x -= split_pos;
        } else if (axis == 1) {
            c1.size.y = split_pos;
            c2.min.y += split_pos;
            c2.size.y -= split_pos;
        } else {
            c1.size.z = split_pos;
            c2.min.z += split_pos;
            c2.size.z -= split_pos;
        }
        
        // Replace parent with children
        leaves.erase(leaves.begin() + best_idx);
        leaves.push_back(c1);
        leaves.push_back(c2);
    }
    
    // 3. Assign Rooms (Shuffle leaves for randomness)
    // We have `leaves` (potentially more than rooms, or eq).
    // Assign processing_rooms[i] to leaves[i].
    
    // Helper to shuffle:
    // ... rng based shuffle manually ...
    for (int i=0; i<leaves.size(); ++i) {
        int swap_idx = rng->randi_range(0, leaves.size() - 1);
        GridAABB temp = leaves[i];
        leaves[i] = leaves[swap_idx];
        leaves[swap_idx] = temp;
    }

    for (int i = 0; i < processing_rooms.size(); ++i) {
        if (i >= leaves.size()) {
            UtilityFunctions::printerr("BSP: Not enough space/leaves for room ", i);
            break;
        }
        
        GridAABB leaf = leaves[i];
        ResolvedRoom &room = processing_rooms[i];
        
        // fixed constraints override BSP? 
        if (room.is_fixed) {
             // Respect fixed
             Vector3 fix_center = room.fixed_position;
             room.position = Vector3i(
                fix_center.x - room.size.x / 2,
                fix_center.y - room.size.y / 2,
                fix_center.z - room.size.z / 2
            );
             // Should we clamp fixed to bounds? Original code did?
             // Assuming fixed is fixed.
             continue;
        }
        
        // Place room in leaf
        // Room Size vs Leaf Size
        // If room fits, random place. If not, center and clamp (might clip neighboring leaf/margin?).
        // Strictly enforcing margin: Clamp to Leaf Bounds.
        
        int free_x = leaf.size.x - room.size.x;
        int free_y = leaf.size.y - room.size.y;
        int free_z = leaf.size.z - room.size.z;
        
        int offset_x = (free_x > 0) ? rng->randi_range(0, free_x) : free_x / 2;
        int offset_y = (free_y > 0) ? rng->randi_range(0, free_y) : free_y / 2;
        int offset_z = (free_z > 0) ? rng->randi_range(0, free_z) : free_z / 2;
        
        room.position = leaf.min + Vector3i(offset_x, offset_y, offset_z);
        
        // Ensure strictly within margin?
        // Leaf Min is >= 1.
        // Leaf Max (Min+Size) is <= gsx-1.
        // Room Pos = Leaf Min + Offset.
        // Room Max = Room Pos + Size.
        // If offset_x positive (room fits), Room Max = Leaf Min + Offset + Size <= Leaf Min + (Leaf Size - Size) + Size = Leaf Min + Leaf Size = Leaf Max. OK.
        // If offset_x negative (room too big), we deliberately center it, potentially overlapping?
        // "One voxel margin" is hard requirement?
        // Let's strict clamp to World Margin [1, gsx-2]
        
        room.position.x = Math::clamp(room.position.x, 1, gsx - 1 - room.size.x);
        room.position.y = Math::clamp(room.position.y, 1, gsy - 1 - room.size.y);
        room.position.z = Math::clamp(room.position.z, 1, gsz - 1 - room.size.z);
    }
}


Dictionary LevelDensityGrid::get_surface_zones() const {
    std::map<int, std::vector<Vector3>> temp_map;
    Array keys = surface_normals.keys();
    for(int i=0; i<keys.size(); ++i) {
        Vector3i pos = keys[i];
        int z_id = get_zone_at(pos);
        temp_map[z_id].push_back(Vector3(pos));
    }

    Dictionary result;
    for(auto const& [id, vec] : temp_map) {
        PackedVector3Array pva;
        pva.resize(vec.size());
        // memcpy might be faster but vectors are safe
        for(size_t k=0; k<vec.size(); ++k) {
            pva[k] = vec[k];
        }
        result[id] = pva;
    }
    return result;
}


void LevelDensityGrid::populate_multimesh(Object* p_multimesh, const PackedVector3Array& points, float scale) {
    if (!p_multimesh) {
        UtilityFunctions::printerr("LevelDensityGrid::populate_multimesh: passed object is null.");
        return;
    }
    MultiMesh* mm = Object::cast_to<MultiMesh>(p_multimesh);
    if (!mm) {
        UtilityFunctions::printerr("LevelDensityGrid::populate_multimesh: passed object is not a MultiMesh. It is: ", p_multimesh->get_class());
        return;
    }
    
    int count = points.size();
    mm->set_instance_count(count);
    
    for(int i=0; i<count; ++i) {
        Vector3 pos = (points[i] + Vector3(0.5f, 0.5f, 0.5f)) * scale;
        Transform3D t;
        t.origin = pos;
        mm->set_instance_transform(i, t);
    }
}

} // namespace godot