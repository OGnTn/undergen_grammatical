#include "path_carver.h"
#include <godot_cpp/variant/utility_functions.hpp>
#include <godot_cpp/core/math.hpp>
#include <godot_cpp/classes/time.hpp>

namespace godot {

void PathCarver::connect_rooms(DensityGrid* grid, RandomNumberGenerator* rng, Ref<FastNoiseLite> wobble_noise, const TypedArray<Dictionary>& generated_rooms, float voxel_size, Vector3& out_spawn, Vector3& out_end) {
    if (generated_rooms.size() < 2) {
        UtilityFunctions::printerr("Not enough rooms generated (", generated_rooms.size(), ") to create paths.");
        if (grid->get_grid_size_x() > 0) {
            Vector3i center = grid->get_grid_dimensions() / 2;
            Vector3i spawn_cell = _find_ground_position(grid, center);
            out_spawn = (Vector3(spawn_cell) + Vector3(0.5f, 0.5f, 0.5f)) * voxel_size;
            out_end = out_spawn;
        } else {
            out_spawn = Vector3(0,0,0);
            out_end = Vector3(0,0,0);
        }
        return;
    }
    
    Time* time = Time::get_singleton();
    uint64_t start_time = time->get_ticks_usec();

    if (dungeon_mode) {
        UtilityFunctions::print("Connecting ", generated_rooms.size(), " rooms with dungeon-style rectangular corridors.");
    } else {
        UtilityFunctions::print("Connecting ", generated_rooms.size(), " rooms with bezier curve paths.");
    }

    const float max_x = (float)grid->get_grid_size_x() - 1.0f;
    const float max_y = (float)grid->get_grid_size_y() - 1.0f;
    const float max_z = (float)grid->get_grid_size_z() - 1.0f;

    for (int i = 0; i < generated_rooms.size() - 1; ++i) {
        Dictionary start_room = generated_rooms[i];
        Dictionary end_room = generated_rooms[i + 1];

        Vector3i start_room_start = start_room["start"];
        Vector3i start_room_end = start_room["end"];
        Vector3i end_room_start = end_room["start"];
        Vector3i end_room_end = end_room["end"];

        Vector3 path_start_point;
        Vector3 path_end_point;

        if (connect_from_ground_level) {
            Vector3 start_room_center = Vector3(start_room_start + (start_room_end - start_room_start) / 2);
            Vector3 end_room_center = Vector3(end_room_start + (end_room_end - end_room_start) / 2);
            Vector3 direction_vector = end_room_center - start_room_center;

            int start_y = Math::max(start_room_start.y, 1);
            int start_x, start_z;
            if (abs(direction_vector.x) > abs(direction_vector.z)) {
                start_z = rng->randi_range(start_room_start.z, start_room_end.z - 1);
                start_x = (direction_vector.x > 0) ? start_room_end.x - 1 : start_room_start.x;
            } else {
                start_x = rng->randi_range(start_room_start.x, start_room_end.x - 1);
                start_z = (direction_vector.z > 0) ? start_room_end.z - 1 : start_room_start.z;
            }
            path_start_point = Vector3(start_x, start_y, start_z);

            int end_y = Math::max(end_room_start.y, 1);
            int end_x, end_z;
             if (abs(direction_vector.x) > abs(direction_vector.z)) {
                end_z = rng->randi_range(end_room_start.z, end_room_end.z - 1);
                end_x = (direction_vector.x > 0) ? end_room_start.x : end_room_end.x - 1;
            } else {
                end_x = rng->randi_range(end_room_start.x, end_room_end.x - 1);
                end_z = (direction_vector.z > 0) ? end_room_start.z : end_room_end.z - 1;
            }
            path_end_point = Vector3(end_x, end_y, end_z);

        } else {
            path_start_point = Vector3(start_room_start + (start_room_end - start_room_start) / 2);
            path_end_point = Vector3(end_room_start + (end_room_end - end_room_start) / 2);
        }

        if (dungeon_mode) {
            if (dungeon_path_algorithm == ALGO_CASTLE_RECURSIVE) {
                _carve_path_castle(grid, rng, path_start_point, path_end_point);
            } else {
                _carve_dungeon_path(grid, path_start_point, path_end_point);
            }
        } else {
            // BEZIER LOGIC
            // NOTE: _mark_brush requires grid and rng.
            std::vector<Vector3> path_nodes;
            path_nodes.push_back(path_start_point);
            int segments = (path_segments > 1) ? path_segments : 1;
            if (segments > 1) {
                Vector3 overall_delta = path_end_point - path_start_point;
                for (int j = 1; j < segments; ++j) {
                    float fraction = static_cast<float>(j) / static_cast<float>(segments);
                    Vector3 base_point = path_start_point + overall_delta * fraction;
                    float offset_mag = (overall_delta.length() / static_cast<float>(segments)) * 0.75f;
                    Vector3 offset(rng->randf_range(-offset_mag, offset_mag), rng->randf_range(-offset_mag, offset_mag), rng->randf_range(-offset_mag, offset_mag));
                    Vector3 intermediate_node = base_point + offset;
                    intermediate_node.x = Math::clamp(intermediate_node.x, 0.0f, max_x);
                    intermediate_node.y = Math::clamp(intermediate_node.y, 0.0f, max_y);
                    intermediate_node.z = Math::clamp(intermediate_node.z, 0.0f, max_z);
                    path_nodes.push_back(intermediate_node);
                }
            }
            path_nodes.push_back(path_end_point);

            for (size_t j = 0; j < path_nodes.size() - 1; ++j) {
                Vector3 p0 = path_nodes[j];
                Vector3 p3 = path_nodes[j + 1];
                Vector3 delta = p3 - p0;
                float offset_magnitude = delta.length() * path_bend_factor;
                Vector3 cp1_offset(rng->randf_range(-offset_magnitude, offset_magnitude), rng->randf_range(-offset_magnitude, offset_magnitude), rng->randf_range(-offset_magnitude, offset_magnitude));
                Vector3 cp2_offset(rng->randf_range(-offset_magnitude, offset_magnitude), rng->randf_range(-offset_magnitude, offset_magnitude), rng->randf_range(-offset_magnitude, offset_magnitude));
                Vector3 p1 = p0 + (delta / 3.0f) + cp1_offset;
                Vector3 p2 = p0 + (delta * 2.0f / 3.0f) + cp2_offset;
                p1.x = Math::clamp(p1.x, 0.0f, max_x); p1.y = Math::clamp(p1.y, 0.0f, max_y); p1.z = Math::clamp(p1.z, 0.0f, max_z);
                p2.x = Math::clamp(p2.x, 0.0f, max_x); p2.y = Math::clamp(p2.y, 0.0f, max_y); p2.z = Math::clamp(p2.z, 0.0f, max_z);
                int num_steps = static_cast<int>(delta.length());
                if (num_steps < 10) num_steps = 10;
                for (int k = 0; k <= num_steps; ++k) {
                    float t = static_cast<float>(k) / static_cast<float>(num_steps);
                    float one_minus_t = 1.0f - t; float t2 = t * t; float t3 = t2 * t; float omt2 = one_minus_t * one_minus_t; float omt3 = omt2 * one_minus_t;
                    Vector3 point_on_curve = (p0 * omt3) + (p1 * 3.0f * omt2 * t) + (p2 * 3.0f * one_minus_t * t2) + (p3 * t3);
                    Vector3 final_point = point_on_curve;
                    
                    if (wobble_noise.is_valid() && path_wobble_magnitude > 0.001f) {
                        final_point.x += wobble_noise->get_noise_1d(point_on_curve.x * path_wobble_frequency) * path_wobble_magnitude;
                        final_point.y += wobble_noise->get_noise_1d(point_on_curve.y * path_wobble_frequency + 1000.0f) * path_wobble_magnitude;
                        final_point.z += wobble_noise->get_noise_1d(point_on_curve.z * path_wobble_frequency + 2000.0f) * path_wobble_magnitude;
                    }
                    _mark_brush(grid, rng, Vector3i(final_point), path_brush_min_radius, path_brush_max_radius, WORLD_OPEN_VALUE);
                }
            }
        }
    }

    // Calc spawn endpoints
    if (generated_rooms.size() > 0) {
        Dictionary first_room = generated_rooms[0];
        Dictionary last_room = generated_rooms[generated_rooms.size() - 1];
        Vector3i start_room_start = first_room["start"];
        Vector3i start_room_end = first_room["end"];
        Vector3i start_spawn_midpoint = start_room_start + (start_room_end - start_room_start) / 2;
        Vector3i spawn_cell = _find_ground_position(grid, start_spawn_midpoint);
        out_spawn = (Vector3(spawn_cell) + Vector3(0.5f, 0.5f, 0.5f)) * voxel_size;
        
        Vector3i end_room_start = last_room["start"];
        Vector3i end_room_end = last_room["end"];
        Vector3i end_room_midpoint = end_room_start + (end_room_end - end_room_start) / 2;
        Vector3i end_cell = _find_ground_position(grid, end_room_midpoint);
        out_end = (Vector3(end_cell) + Vector3(0.5f, 0.5f, 0.5f)) * voxel_size;
    }
    
    uint64_t end_time = time->get_ticks_usec();
    double total_ms = (end_time - start_time) / 1000.0;
    UtilityFunctions::print("PathCarver::connect_rooms Total Time: ", total_ms, " ms");
}

// Updated signature
void PathCarver::create_paths_from_edges(DensityGrid* grid, RandomNumberGenerator* rng, Ref<FastNoiseLite> wobble_noise, const std::vector<ResolvedRoom>& rooms, const std::vector<ResolvedEdge>& edges) {
    for(const auto &edge : edges) {
        int z_id = grid->register_zone_name(edge.type);
        current_carving_zone_id = z_id;

        // Safety Check for Room Indices
        if (edge.from_index < 0 || edge.from_index >= rooms.size() || 
            edge.to_index < 0 || edge.to_index >= rooms.size()) {
            UtilityFunctions::printerr("PathCarver: Edge references invalid room indices: ", edge.from_index, " -> ", edge.to_index, ". Max: ", rooms.size());
            continue;
        }

        const ResolvedRoom &rA = rooms[edge.from_index];
        const ResolvedRoom &rB = rooms[edge.to_index];
        
        Vector3 start_point;
        Vector3 end_point;

        if (connect_from_ground_level) {
            Vector3i rA_start = rA.position;
            Vector3i rA_end = rA.position + rA.size;
            Vector3i rB_start = rB.position;
            Vector3i rB_end = rB.position + rB.size;

            Vector3 rA_center = rA.center();
            Vector3 rB_center = rB.center();
            Vector3 direction_vector = rB_center - rA_center;

            // Start Point Calculation
            int start_y = Math::max(rA_start.y, 1);
            int start_x, start_z;
            
            if (abs(direction_vector.x) > abs(direction_vector.z)) {
                // Moving primarily along X
                start_z = rng->randi_range(rA_start.z, rA_end.z - 1);
                start_x = (direction_vector.x > 0) ? rA_end.x - 1 : rA_start.x;
            } else {
                // Moving primarily along Z
                start_x = rng->randi_range(rA_start.x, rA_end.x - 1);
                start_z = (direction_vector.z > 0) ? rA_end.z - 1 : rA_start.z;
            }
            start_point = Vector3(start_x, start_y, start_z);

            // End Point Calculation
            int end_y = Math::max(rB_start.y, 1);
            int end_x, end_z;
            
            if (abs(direction_vector.x) > abs(direction_vector.z)) {
                // Moving primarily along X
                end_z = rng->randi_range(rB_start.z, rB_end.z - 1);
                // If moving +X, we enter from Left (start_x). If moving -X, we enter from Right (end_x-1)
                end_x = (direction_vector.x > 0) ? rB_start.x : rB_end.x - 1;
            } else {
                // Moving primarily along Z
                end_x = rng->randi_range(rB_start.x, rB_end.x - 1);
                end_z = (direction_vector.z > 0) ? rB_start.z : rB_end.z - 1;
            }
            end_point = Vector3(end_x, end_y, end_z);

        } else {
            start_point = rA.center();
            end_point = rB.center();
        }

        if(dungeon_mode) {
             _carve_dungeon_path(grid, start_point, end_point);
        } else {
             _carve_bezier_path(grid, rng, wobble_noise, start_point, end_point);
        }
    }
    current_carving_zone_id = 0;
}

void PathCarver::_carve_bezier_path(DensityGrid* grid, RandomNumberGenerator* rng, Ref<FastNoiseLite> wobble_noise, const Vector3 &start, const Vector3 &end) {
    std::vector<Vector3> path_nodes;
    path_nodes.push_back(start);
    
    int segments = (path_segments > 1) ? path_segments : 1;
    float max_x = (float)grid->get_grid_size_x() - 1.0f;
    float max_y = (float)grid->get_grid_size_y() - 1.0f;
    float max_z = (float)grid->get_grid_size_z() - 1.0f;

    if (segments > 1) {
        Vector3 overall_delta = end - start;
        for (int j = 1; j < segments; ++j) {
            float fraction = static_cast<float>(j) / static_cast<float>(segments);
            Vector3 base_point = start + overall_delta * fraction;
            
            // Random Offset for natural look
            float offset_mag = (overall_delta.length() / static_cast<float>(segments)) * 0.75f;
            Vector3 offset(
                rng->randf_range(-offset_mag, offset_mag), 
                rng->randf_range(-offset_mag, offset_mag), 
                rng->randf_range(-offset_mag, offset_mag)
            );
            
            Vector3 intermediate_node = base_point + offset;
            intermediate_node.x = Math::clamp(intermediate_node.x, 0.0f, max_x);
            intermediate_node.y = Math::clamp(intermediate_node.y, 0.0f, max_y);
            intermediate_node.z = Math::clamp(intermediate_node.z, 0.0f, max_z);
            path_nodes.push_back(intermediate_node);
        }
    }
    path_nodes.push_back(end);

    for (size_t j = 0; j < path_nodes.size() - 1; ++j) {
        Vector3 p0 = path_nodes[j];
        Vector3 p3 = path_nodes[j + 1];
        Vector3 delta = p3 - p0;
        
        float offset_magnitude = delta.length() * path_bend_factor;
        
        Vector3 cp1_offset(
            rng->randf_range(-offset_magnitude, offset_magnitude), 
            rng->randf_range(-offset_magnitude, offset_magnitude), 
            rng->randf_range(-offset_magnitude, offset_magnitude)
        );
        Vector3 cp2_offset(
            rng->randf_range(-offset_magnitude, offset_magnitude), 
            rng->randf_range(-offset_magnitude, offset_magnitude), 
            rng->randf_range(-offset_magnitude, offset_magnitude)
        );
        
        Vector3 p1 = p0 + (delta / 3.0f) + cp1_offset;
        Vector3 p2 = p0 + (delta * 2.0f / 3.0f) + cp2_offset;
        
        // Clamp Control Points
        p1.x = Math::clamp(p1.x, 0.0f, max_x); p1.y = Math::clamp(p1.y, 0.0f, max_y); p1.z = Math::clamp(p1.z, 0.0f, max_z);
        p2.x = Math::clamp(p2.x, 0.0f, max_x); p2.y = Math::clamp(p2.y, 0.0f, max_y); p2.z = Math::clamp(p2.z, 0.0f, max_z);
        
        int num_steps = static_cast<int>(delta.length());
        if (num_steps < 10) num_steps = 10;
        
        for (int k = 0; k <= num_steps; ++k) {
            float t = static_cast<float>(k) / static_cast<float>(num_steps);
            float one_minus_t = 1.0f - t; 
            float t2 = t * t; 
            float t3 = t2 * t; 
            float omt2 = one_minus_t * one_minus_t; 
            float omt3 = omt2 * one_minus_t;
            
            // Cubic Bezier validation
            Vector3 point_on_curve = (p0 * omt3) + (p1 * 3.0f * omt2 * t) + (p2 * 3.0f * one_minus_t * t2) + (p3 * t3);
            Vector3 final_point = point_on_curve;
            
            if (path_wobble_magnitude > 0.001f && wobble_noise.is_valid()) {
                final_point.x += wobble_noise->get_noise_1d(point_on_curve.x * path_wobble_frequency) * path_wobble_magnitude;
                final_point.y += wobble_noise->get_noise_1d(point_on_curve.y * path_wobble_frequency + 1000.0f) * path_wobble_magnitude;
                final_point.z += wobble_noise->get_noise_1d(point_on_curve.z * path_wobble_frequency + 2000.0f) * path_wobble_magnitude;
            }
            
            // Varying width noise
            int eff_radius = path_brush_min_radius;
            if (cave_width_noise > 0.001f && wobble_noise.is_valid()) {
                 float w_noise = wobble_noise->get_noise_1d(static_cast<float>(k) * 0.1f + j * 10.0f); // vary along path
                 eff_radius = (int)((float)path_brush_min_radius * (1.0f + w_noise * cave_width_noise));
                 if (eff_radius < 1) eff_radius = 1;
            } else {
                 // Use range if no varying noise, otherwise single radius
                 eff_radius = (path_brush_min_radius == path_brush_max_radius) ? path_brush_min_radius : rng->randi_range(path_brush_min_radius, path_brush_max_radius);
            }

            _carve_complex_brush(grid, Vector3i(final_point), eff_radius, wobble_noise);
        }
    }
}

void PathCarver::_carve_complex_brush(DensityGrid* grid, const Vector3i &center, int radius, Ref<FastNoiseLite> noise) {
    if (radius <= 0) return;
    
    // Bounds check optimization
    int gsx = grid->get_grid_size_x();
    int gsy = grid->get_grid_size_y();
    int gsz = grid->get_grid_size_z();
    
    if (center.x + radius < 0 || center.x - radius >= gsx ||
        center.y + radius < 0 || center.y - radius >= gsy ||
        center.z + radius < 0 || center.z - radius >= gsz) {
        return;
    }

    // Iterate bounding box
    // To handle "ruggedness" effectively, we might need a slightly larger bounding box than just 'radius' 
    // if the noise adds to it. But generally 'radius' is the target average.
    // Let's pad slightly if ruggedness is high.
    int pad = (int)Math::ceil(Math::max(cave_ruggedness, Math::max(cave_floor_ruggedness, cave_ceiling_ruggedness)) * 2.0f);
    int iter_rad = radius + pad;

    for (int i = -iter_rad; i <= iter_rad; ++i) {
        for (int j = -iter_rad; j <= iter_rad; ++j) {
            for (int k = -iter_rad; k <= iter_rad; ++k) {
                Vector3 offset(i, j, k);
                Vector3i pos = center + Vector3i(i, j, k);

                if (!grid->is_valid_position(pos)) continue;

                // 1. Base Distance
                float dist = offset.length();

                // 2. Overhang (Modify apparent vertical distance)
                // If y > 0, we shrink the apparent vertical capability to make it look wider?
                // Actually, to make it wider at top, we want 'dist' to be smaller for same offset.
                // Keyhole: Top is wide, Bottom is narrow.
                // If j > 0 (top), we multiply j by something < 1.0 (squash) -> implies wider shape
                // OR we can just modify the target radius.
                float current_radius_target = (float)radius;
                
                if (overhang_openness > 0.001f) {
                   // Trapezoidal logic: 
                   // If j > 0, radius expands.
                   // If j < 0, radius shrinks?
                   float rel_y = (float)j / (float)radius; // -1 to 1
                   // curve: 1.0 + rel_y * factor
                   float shape_factor = 1.0f + (rel_y * overhang_openness); 
                   // Clamp to avoid inverted shapes if extreme
                   if (shape_factor < 0.1f) shape_factor = 0.1f;
                   
                   current_radius_target *= shape_factor;
                }

                // 3. Noise / Ruggedness
                if (noise.is_valid()) {
                    // Use world coordinates for noise continuity
                    float n_val = noise->get_noise_3d(pos.x, pos.y, pos.z); 
                    
                    // General wall ruggedness
                    current_radius_target += n_val * cave_ruggedness * 3.0f; // * 3.0 arbitrary scaling
                    
                    // Specific Floor/Ceiling
                    if (j < -radius/3) { // Floor area
                         float n_floor = noise->get_noise_3d(pos.x, pos.y * 2.0f, pos.z);
                         // Stalagmites: noise reduces radius (grows into cave)
                         // If n_floor is positive -> adds to radius (hole gets bigger)
                         // If n_floor is negative -> subtracts (rock stays)
                         // We want stalagmites (rock sticking up). So we want to SUBTRACT from radius.
                         // But noise is -1..1.
                         current_radius_target -= Math::abs(n_floor) * cave_floor_ruggedness * 4.0f;
                    } 
                    else if (j > radius/3) { // Ceiling area
                         float n_ceil = noise->get_noise_3d(pos.x, pos.y * 2.0f + 500, pos.z);
                         // Stalactites
                         current_radius_target -= Math::abs(n_ceil) * cave_ceiling_ruggedness * 4.0f;
                    }
                }

                // 4. Floor Flattening
                // If we are deep enough (low Y), we forcefully mask out.
                if (floor_flattening > 0.001f) {
                    float bottom_threshold = -(float)radius * (1.0f - floor_flattening);
                    if ((float)j < bottom_threshold) {
                         continue; // Don't carve this (leave as separate ground or solid)
                    }
                }

                // Final check
                if (dist < current_radius_target) {
                    grid->set_cell(pos, WORLD_OPEN_VALUE);
                     // Manual boundary check for Zone safety
                    if (current_carving_zone_id > 0) {
                        if (grid->get_zone_at(pos) == 0) {
                            grid->set_zone_at(pos, current_carving_zone_id);
                        }
                    }
                }
            }
        }
    }
}

// ================= PRIVATE HELPERS =================

void PathCarver::_mark_brush(DensityGrid* grid, RandomNumberGenerator* rng, const Vector3i &center, int radius_lower, int radius_higher, float value) {
    int r_low = radius_lower;
    int r_high = radius_higher;
    if (r_high < r_low) r_high = r_low;
    if (r_low < 0) r_low = 0;

    int radius = (r_low == r_high) ? r_low : rng->randi_range(r_low, r_high);
    if (radius <= 0) return;

    int radius_sq = radius * radius;
    int gsx = grid->get_grid_size_x();
    int gsy = grid->get_grid_size_y();
    int gsz = grid->get_grid_size_z();

    if(dungeon_mode) {
        for (int i = -radius; i <= radius; ++i) {
            for (int j = 0; j <= radius * 2; ++j) { 
                for (int k = -radius; k <= radius; ++k) {
                    // Box Brush
                    Vector3i current_pos = center + Vector3i(i, j, k);
                    
                    // explicit boundary check (margin 1)
                    if (current_pos.x > 0 && current_pos.x < gsx - 1 &&
                        current_pos.y > 0 && current_pos.y < gsy - 1 &&
                        current_pos.z > 0 && current_pos.z < gsz - 1) {
                        
                        grid->set_cell(current_pos, value);
                            // Manual boundary check for Zone safety
                        if (current_carving_zone_id > 0 && grid->is_valid_position(current_pos)) {
                            if (grid->get_zone_at(current_pos) == 0) {
                                grid->set_zone_at(current_pos, current_carving_zone_id);
                            }
                        }
                    }
                }
            }
        }
    }
    else {
        for (int i = -radius; i <= radius; ++i) {
            for (int j = -radius; j <= radius; ++j) {
                for (int k = -radius; k <= radius; ++k) {
                    bool in_shape = false;
                    if (use_square_brush) {
                        in_shape = true; 
                    } else {
                        if (i * i + j * j + k * k < radius_sq) {
                            in_shape = true; 
                        }
                    }

                    if (in_shape) {
                        Vector3i current_pos = center + Vector3i(i, j, k);
                        // explicit boundary check (margin 1)
                        if (current_pos.x > 0 && current_pos.x < gsx - 1 &&
                            current_pos.y > 0 && current_pos.y < gsy - 1 &&
                            current_pos.z > 0 && current_pos.z < gsz - 1) {
                            
                            grid->set_cell(current_pos, value);
                            
                            // Manual boundary check for Zone safety
                            if (current_carving_zone_id > 0 && grid->is_valid_position(current_pos)) {
                                if (grid->get_zone_at(current_pos) == 0) {
                                    grid->set_zone_at(current_pos, current_carving_zone_id);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

Vector3i PathCarver::_find_ground_position(DensityGrid* grid, const Vector3i &start_pos) {
    int gsx = grid->get_grid_size_x();
    int gsy = grid->get_grid_size_y();
    int gsz = grid->get_grid_size_z();
    
    int x = Math::clamp(start_pos.x, 0, gsx - 1);
    int z = Math::clamp(start_pos.z, 0, gsz - 1);
    int start_y = Math::clamp(start_pos.y, 0, gsy - 1);

    // Search Down
    for (int y = start_y; y >= 0; --y) {
        Vector3i pos(x, y, z);
        if (grid->get_cell(pos, WORLD_SOLID_VALUE) >= grid->get_surface_threshold()) {
            return Vector3i(x, y + 1, z);
        }
    }
    return Vector3i(x, 0, z); // Fallback to bottom
}

void PathCarver::_carve_corridor_segment(DensityGrid* grid, RandomNumberGenerator* rng, const Vector3i &from, const Vector3i &to) {
    if (from == to) {
        _mark_brush(grid, rng, from, path_brush_min_radius, path_brush_max_radius, WORLD_OPEN_VALUE);
        return;
    }

    // Manhattan Carving (X -> Z -> Y)
    // This guarantees reaching the target without loops
    Vector3i cursor = from;
    
    // 1. Move X
    int x_dir = (to.x > cursor.x) ? 1 : -1;
    while (cursor.x != to.x) {
        cursor.x += x_dir;
        _mark_brush(grid, rng, cursor, path_brush_min_radius, path_brush_max_radius, WORLD_OPEN_VALUE);
    }
    
    // 2. Move Z
    int z_dir = (to.z > cursor.z) ? 1 : -1;
    while (cursor.z != to.z) {
        cursor.z += z_dir;
        _mark_brush(grid, rng, cursor, path_brush_min_radius, path_brush_max_radius, WORLD_OPEN_VALUE);
    }

    // 3. Move Y
    int y_dir = (to.y > cursor.y) ? 1 : -1;
    while (cursor.y != to.y) {
        cursor.y += y_dir;
        _mark_brush(grid, rng, cursor, path_brush_min_radius, path_brush_max_radius, WORLD_OPEN_VALUE);
    }
}

// ================= DUNGEON A* =================

void PathCarver::_carve_dungeon_path(DensityGrid* grid, const Vector3 &start, const Vector3 &end) {
    Vector3i start_i(start);
    Vector3i end_i(end);

    int gsx = grid->get_grid_size_x();
    int gsy = grid->get_grid_size_y();
    int gsz = grid->get_grid_size_z();

    Time* time = Time::get_singleton();
    uint64_t t_start = time->get_ticks_usec();

    float surf_thresh = grid->get_surface_threshold();

    // Config Constants
    const int COST_MOVE_BASE = 10;
    const int COST_ROOM_PENALTY = 100; // Prefer existing paths
    const int COST_TURN_BASE = 50;
    const int COST_TURN_MIN = 0;
    const float COST_TURN_DECAY = 2.0f; // Reduction per straight step
    
    // Stair Config
    const int COST_STAIR_START_BASE = 80;
    const int COST_STAIR_START_MIN = 10;
    const float COST_STAIR_START_DECAY = 5.0f; 
    const int COST_STAIR_STEP_BASE = (int)(COST_MOVE_BASE * vertical_movement_cost_multiplier);
    const int COST_STAIR_LEN_PENALTY = 2; 

    // Heuristics
    const int HEURISTIC_WEIGHT = 2;

    std::priority_queue<DungeonNode, std::vector<DungeonNode>, std::greater<DungeonNode>> open_set;
    std::unordered_map<int64_t, int> closed_set_g;
    std::unordered_map<int64_t, DungeonNode> came_from;

    auto get_idx = [&](const Vector3i& p) -> int64_t {
        return (int64_t)p.z * gsy * gsx + (int64_t)p.y * gsx + p.x;
    };

    DungeonNode start_node = { start_i, Vector3i(0,0,0), 0, 0, 0, 0, 0, false, start_i };
    open_set.push(start_node);
    closed_set_g[get_idx(start_i)] = 0;

    DungeonNode final_node;
    bool found = false;
    int iterations = 0;
    int max_iterations = 1000000;

    while (!open_set.empty()) {
        DungeonNode current = open_set.top();
        open_set.pop();

        if (iterations > max_iterations) break;
        if (current.pos == end_i) {
            final_node = current;
            found = true;
            break;
        }

        int current_g = closed_set_g[get_idx(current.pos)];
        if (current.g > current_g) continue; 

        // Directions: 4 Cardinal + 2 Vertical
        struct DirOpt { Vector3i dir; bool is_stair; int y_change; };
        
        Vector3i d_flat = current.arrived_from_dir; 
        if (d_flat.y != 0) d_flat = Vector3i(0,0,0); // If we came from vert, we have no flat momentum

        // Available moves
        std::vector<DirOpt> options = {
            { Vector3i(1, 0, 0), false, 0 }, { Vector3i(-1, 0, 0), false, 0 },
            { Vector3i(0, 0, 1), false, 0 }, { Vector3i(0, 0, -1), false, 0 }
        };

        // Add stairs (must move 1 flat + 1 vert)
        // We only add stairs if we have a valid flat direction (momentum) OR we pick a new one?
        // Actually, "Stair" implies moving forward AND up/down.
        // For simplicity from original code: we allow moving in current flat dir + up/down.
        // OR if stationary, try all 4. 
        // Original code hardcoded some specific logic. I will implement simplified "Moves":
        // 1. Move Flat (Cardinal)
        // 2. Move Stair (Cardinal + Up/Down)
        
        // Let's expand options to include stairs for all cardinals
        std::vector<DirOpt> expanded_options;
        for(auto& o : options) {
            expanded_options.push_back(o); // Flat
            expanded_options.push_back({ o.dir + Vector3i(0, 1, 0), true, 1 }); // Up
            expanded_options.push_back({ o.dir + Vector3i(0, -1, 0), true, -1 }); // Down
        }

        for (const auto& opt : expanded_options) {
            Vector3i next_pos = current.pos + opt.dir; // e.g. (1, 1, 0)
            
            // 1. Boundary Check (Relaxed to allow rooms near edge)
            if (next_pos.x < 1 || next_pos.x >= gsx - 1 ||
                next_pos.y < 1 || next_pos.y >= gsy - 1 ||
                next_pos.z < 1 || next_pos.z >= gsz - 1) {
                continue;
            }

            // 2. Determine State
            // Is this a straight move?
            // "Straight" means the flat component matches current.arrived_from_dir
            Vector3i flat_move_dir = Vector3i(opt.dir.x, 0, opt.dir.z); // Extract flat part
            // Note: opt.dir for stair is (1,1,0). Flat part (1,0,0).
            bool is_straight = (flat_move_dir == d_flat) || (d_flat == Vector3i(0,0,0));

            int move_cost = COST_MOVE_BASE;

            // A. Room Avoidance
            float cell_val = grid->get_cell(next_pos, WORLD_SOLID_VALUE);
            if (cell_val < surf_thresh && next_pos != end_i && next_pos != start_i) {
                move_cost += COST_ROOM_PENALTY;
            }

            // B. Turn Logic
            if (!is_straight) {
                int reduction = current.straight_dist * COST_TURN_DECAY;
                int turn_cost = Math::max(COST_TURN_MIN, COST_TURN_BASE - reduction);
                move_cost += turn_cost;
            }

            // C. Stair Logic
             if (opt.is_stair) {
                if (current.is_in_staircase) {
                    // Continuing
                    Vector3i prev_move = current.pos - current.parent_pos;
                    if (Math::sign(prev_move.y) != opt.y_change) {
                        move_cost += 2000; // Penalty for zigzag
                    }
                    move_cost += COST_STAIR_STEP_BASE;
                    move_cost += (current.current_stair_len * COST_STAIR_LEN_PENALTY);
                } else {
                    // Starting
                    int reduction = current.dist_since_stair * COST_STAIR_START_DECAY;
                    int start_cost = Math::max(COST_STAIR_START_MIN, COST_STAIR_START_BASE - reduction);
                    move_cost += start_cost;
                }
            }

            int new_straight_dist = is_straight ? current.straight_dist + 1 : 1;
            int new_dist_since_stair = current.dist_since_stair;
            int new_stair_len = 0;

            if (opt.is_stair) {
                new_stair_len = current.current_stair_len + 1;
                new_dist_since_stair = 0;
            } else {
                new_stair_len = 0;
                if (!current.is_in_staircase) {
                     new_dist_since_stair++;
                }
            }
            
            Vector3i new_entered_from = flat_move_dir;
            
            int new_g = current.g + move_cost;
            int64_t n_idx = get_idx(next_pos);

            if (closed_set_g.find(n_idx) == closed_set_g.end() || new_g < closed_set_g[n_idx]) {
                closed_set_g[n_idx] = new_g;
                
                int new_h = (Math::abs(next_pos.x - end_i.x) + 
                             Math::abs(next_pos.y - end_i.y) + 
                             Math::abs(next_pos.z - end_i.z)) * COST_MOVE_BASE * HEURISTIC_WEIGHT;

                DungeonNode next_node = {
                    next_pos, new_entered_from, new_g, new_h,
                    new_straight_dist, new_dist_since_stair, new_stair_len, 
                    opt.is_stair, current.pos
                };
                came_from[n_idx] = next_node;
                open_set.push(next_node);
            }
        }
        iterations++;
    }

    if (found) {
        // Reconstruct
        std::vector<Vector3i> path;
        DungeonNode curr = final_node;
        while (curr.pos != start_i) {
            path.push_back(curr.pos);
            int64_t idx = get_idx(curr.pos);
            if (came_from.find(idx) == came_from.end()) break;
            
            // Parent Lookup
            int64_t parent_idx = get_idx(curr.parent_pos);
            if (came_from.find(parent_idx) != came_from.end()) {
                curr = came_from[parent_idx];
            } else if (curr.parent_pos == start_i) {
                break;
            } else {
                break;
            }
        }
        path.push_back(start_i);

        // Dummy RNG for brush (if needed) - we'll just pass nullptr if we used a strict mark method,
        // but our _mark_brush uses RNG for jitter.
        // We need an RNG instance. Since we don't have one passed here easily (wait, we have one in class if we stored it?)
        // Design fix: pass RNG to _carve_dungeon_path.
        // I'll create a temp local RNG if needed or update signature.
        // Updating signature is cleaner. For now using a local one to save time/complexity of update?
        // No, I added Update signature in header: `_carve_dungeon_path(DensityGrid*, ...)`
        // I need to update header to pass RNG.
        // Re-reading `path_carver.h` I just wrote...
        // `_carve_dungeon_path(DensityGrid* grid, const Vector3 &start, const Vector3 &end);`
        // I should have passed RNG. I will instantiate a local one to avoid header rewrite, 
        // OR better: use `RandomNumberGenerator rng; rng.set_seed(1234);`
        // Actually, keeping strict path visualization is better without jitter.
        // So I'll pass a dummy RNG or create one.
        Ref<RandomNumberGenerator> temp_rng;
        temp_rng.instantiate();
        temp_rng->set_seed(12345);

        for (const Vector3i& p : path) {
            _mark_brush(grid, temp_rng.ptr(), p, path_brush_min_radius, path_brush_max_radius, WORLD_OPEN_VALUE);
            grid->set_cell(p + Vector3i(0, 1, 0), WORLD_OPEN_VALUE);
            grid->set_cell(p + Vector3i(0, 2, 0), WORLD_OPEN_VALUE);
        }
    } else {
        UtilityFunctions::printerr("Dungeon Pathfinding Failed. Falling back.");
        Ref<RandomNumberGenerator> temp_rng;
        temp_rng.instantiate();
        _carve_corridor_segment(grid, temp_rng.ptr(), start_i, end_i);
    }
}

// ================= CASTLE / RECURSIVE =================

void PathCarver::_carve_path_castle(DensityGrid* grid, RandomNumberGenerator* rng, const Vector3 &start, const Vector3 &end) {
    _carve_recursive_winding_path(grid, rng, Vector3i(start), Vector3i(end), 4);
}

void PathCarver::_carve_recursive_winding_path(DensityGrid* grid, RandomNumberGenerator* rng, const Vector3i &start, const Vector3i &end, int depth) {
    Vector3 diff = Vector3(end - start);
    float dist = diff.length();

    if (depth <= 0 || dist < 20.0f) {
        _carve_stepped_L_shape(grid, rng, start, end);
        return;
    }

    Vector3i mid = (start + end) / 2;
    mid.y = (start.y + end.y) / 2;

    int offset_mag = (int)(dist * 0.3f); 
    if (offset_mag < 2) offset_mag = 2;

    if (Math::abs(diff.x) > Math::abs(diff.z)) {
        mid.z += rng->randi_range(-offset_mag, offset_mag);
    } else {
        mid.x += rng->randi_range(-offset_mag, offset_mag);
    }
    
    int gsx = grid->get_grid_size_x(); 
    int gsz = grid->get_grid_size_z();
    mid.x = Math::clamp(mid.x, 2, gsx - 3);
    mid.z = Math::clamp(mid.z, 2, gsz - 3);
    
    _carve_recursive_winding_path(grid, rng, start, mid, depth - 1);
    _carve_recursive_winding_path(grid, rng, mid, end, depth - 1);
}

void PathCarver::_carve_stepped_L_shape(DensityGrid* grid, RandomNumberGenerator* rng, const Vector3i &start, const Vector3i &end) {
    bool x_first = rng->randf() > 0.5f;
    bool corner_at_start_height = rng->randf() > 0.5f;

    Vector3i corner_flat;
    if (x_first) {
        corner_flat = Vector3i(end.x, 0, start.z);
    } else {
        corner_flat = Vector3i(start.x, 0, end.z);
    }

    Vector3i corner;
    if (corner_at_start_height) {
        corner = Vector3i(corner_flat.x, start.y, corner_flat.z);
        _carve_variable_height_leg(grid, rng, start, corner);
        _carve_variable_height_leg(grid, rng, corner, end);
    } else {
        corner = Vector3i(corner_flat.x, end.y, corner_flat.z);
        _carve_variable_height_leg(grid, rng, start, corner);
        _carve_variable_height_leg(grid, rng, corner, end);
    }
}

void PathCarver::_carve_variable_height_leg(DensityGrid* grid, RandomNumberGenerator* rng, const Vector3i &start, const Vector3i &end) {
    int height_diff = end.y - start.y;
    
    if (height_diff == 0) {
        _carve_corridor_segment(grid, rng, start, end);
        return;
    }
    
    int h_dist = (int)Vector2(start.x - end.x, start.z - end.z).length();
    int stair_len = Math::abs(height_diff); 
    
    if (stair_len >= h_dist) {
        _carve_corridor_segment(grid, rng, start, end); 
        return;
    }
    
    int flat_padding = (h_dist - stair_len) / 2;
    Vector3i stair_start = start;
    if (start.x != end.x) {
        stair_start.x += Math::sign(end.x - start.x) * flat_padding;
    } else {
        stair_start.z += Math::sign(end.z - start.z) * flat_padding;
    }
    
    Vector3i stair_end = stair_start;
    if (start.x != end.x) {
        stair_end.x += Math::sign(end.x - start.x) * stair_len;
    } else {
        stair_end.z += Math::sign(end.z - start.z) * stair_len;
    }
    stair_end.y = end.y; 

    _carve_corridor_segment(grid, rng, start, stair_start);
    
    // Carve Stairs
    Vector3i cursor = stair_start;
    int y_dir = Math::sign(height_diff);
    Vector3i step_dir = Vector3i(0,0,0);
    if (stair_start.x != stair_end.x) step_dir.x = Math::sign(stair_end.x - stair_start.x);
    else step_dir.z = Math::sign(stair_end.z - stair_start.z);

    int radius = path_brush_min_radius;
    if (radius < 2) radius = 2;

    for (int i = 0; i < stair_len; ++i) {
        _mark_brush(grid, rng, cursor, radius, radius, WORLD_OPEN_VALUE);
        if (radius < 3) _mark_brush(grid, rng, cursor + Vector3i(0, 2, 0), radius, radius, WORLD_OPEN_VALUE);
        cursor.y += y_dir;
        cursor += step_dir;
    }
    _mark_brush(grid, rng, cursor, radius, radius, WORLD_OPEN_VALUE);

    _carve_corridor_segment(grid, rng, stair_end, end);
}

} // namespace godot
