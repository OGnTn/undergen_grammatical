#include "undergen_world_3d.h"
#include "undergen_mesher_node.h"
#include "undergen_scene_spawner_node.h"
#include "undergen_mesh_spawner_node.h"
#include "undergen_grammar_node.h"
#include "undergen_bsp_placer_node.h"
#include "density_grid.h"
#include "undergen_point_set.h"
#include "undergen_vox_stamp_node.h"
#include "undergen_detail_stamper_node.h"
#include "undergen_modular_astar_carver_node.h"
#include "undergen_spatial_model.h"
#include "undergen_spatial_nodes.h"
#include "mc_chunk.h"
#include "dc_chunk.h"
#include <godot_cpp/variant/utility_functions.hpp>
#include <godot_cpp/classes/label3d.hpp>
#include <godot_cpp/classes/multiplayer_spawner.hpp>

namespace godot {

UnderGenWorld3D::UnderGenWorld3D() : is_generating(false) {}

UnderGenWorld3D::~UnderGenWorld3D() {
    cancel_generation();
}

void UnderGenWorld3D::_notification(int p_what) {
    if (p_what == NOTIFICATION_READY) {
        if (generate_on_ready) {
            generate();
        }
    }
}

void UnderGenWorld3D::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_pipeline", "pipeline"), &UnderGenWorld3D::set_pipeline);
    ClassDB::bind_method(D_METHOD("get_pipeline"), &UnderGenWorld3D::get_pipeline);

    ClassDB::bind_method(D_METHOD("set_grammar_override", "grammar"), &UnderGenWorld3D::set_grammar_override);
    ClassDB::bind_method(D_METHOD("get_grammar_override"), &UnderGenWorld3D::get_grammar_override);

    ClassDB::bind_method(D_METHOD("set_generation_seed", "seed"), &UnderGenWorld3D::set_generation_seed);
    ClassDB::bind_method(D_METHOD("get_generation_seed"), &UnderGenWorld3D::get_generation_seed);

    ClassDB::bind_method(D_METHOD("set_voxel_size", "size"), &UnderGenWorld3D::set_voxel_size);
    ClassDB::bind_method(D_METHOD("get_voxel_size"), &UnderGenWorld3D::get_voxel_size);

    ClassDB::bind_method(D_METHOD("set_surface_threshold", "threshold"), &UnderGenWorld3D::set_surface_threshold);
    ClassDB::bind_method(D_METHOD("get_surface_threshold"), &UnderGenWorld3D::get_surface_threshold);

    ClassDB::bind_method(D_METHOD("set_generate_on_ready", "enabled"), &UnderGenWorld3D::set_generate_on_ready);
    ClassDB::bind_method(D_METHOD("get_generate_on_ready"), &UnderGenWorld3D::get_generate_on_ready);

    ClassDB::bind_method(D_METHOD("set_spawn_on_generation_complete", "enabled"), &UnderGenWorld3D::set_spawn_on_generation_complete);
    ClassDB::bind_method(D_METHOD("get_spawn_on_generation_complete"), &UnderGenWorld3D::get_spawn_on_generation_complete);
    ClassDB::bind_method(D_METHOD("set_cast_shadows", "cast_shadows"), &UnderGenWorld3D::set_cast_shadows);
    ClassDB::bind_method(D_METHOD("get_cast_shadows"), &UnderGenWorld3D::get_cast_shadows);
    ClassDB::bind_method(D_METHOD("set_parent_node", "parent_node"), &UnderGenWorld3D::set_parent_node);
    ClassDB::bind_method(D_METHOD("get_parent_node"), &UnderGenWorld3D::get_parent_node);
    ClassDB::bind_method(D_METHOD("set_multiplayer_spawner", "path"), &UnderGenWorld3D::set_multiplayer_spawner);
    ClassDB::bind_method(D_METHOD("get_multiplayer_spawner"), &UnderGenWorld3D::get_multiplayer_spawner);
    ClassDB::bind_method(D_METHOD("spawn_scenes", "parent_node"), &UnderGenWorld3D::spawn_scenes, DEFVAL(nullptr));
    ClassDB::bind_method(D_METHOD("spawn_scenes_for_node", "node_name", "parent_node"), &UnderGenWorld3D::spawn_scenes_for_node, DEFVAL(nullptr));
    ClassDB::bind_method(D_METHOD("get_point_set_from_node", "node_name"), &UnderGenWorld3D::get_point_set_from_node);
    ClassDB::bind_method(D_METHOD("get_vox_spawns"), &UnderGenWorld3D::get_vox_spawns);

    // Inspector button
    ClassDB::bind_method(D_METHOD("set_trigger_generate", "value"), &UnderGenWorld3D::set_trigger_generate);
    ClassDB::bind_method(D_METHOD("get_trigger_generate"), &UnderGenWorld3D::get_trigger_generate);

    // Debug Visualization
    ClassDB::bind_method(D_METHOD("set_debug_show_zone_labels", "enabled"), &UnderGenWorld3D::set_debug_show_zone_labels);
    ClassDB::bind_method(D_METHOD("get_debug_show_zone_labels"), &UnderGenWorld3D::get_debug_show_zone_labels);
    ClassDB::bind_method(D_METHOD("set_debug_zone_label_font_size", "size"), &UnderGenWorld3D::set_debug_zone_label_font_size);
    ClassDB::bind_method(D_METHOD("get_debug_zone_label_font_size"), &UnderGenWorld3D::get_debug_zone_label_font_size);
    ClassDB::bind_method(D_METHOD("set_debug_zone_label_color", "color"), &UnderGenWorld3D::set_debug_zone_label_color);
    ClassDB::bind_method(D_METHOD("get_debug_zone_label_color"), &UnderGenWorld3D::get_debug_zone_label_color);

    // Manual debug label control (callable from GDScript)
    ClassDB::bind_method(D_METHOD("spawn_debug_zone_labels", "context"), &UnderGenWorld3D::_spawn_debug_zone_labels);
    ClassDB::bind_method(D_METHOD("clear_debug_labels"), &UnderGenWorld3D::_clear_debug_labels);

    ClassDB::bind_method(D_METHOD("generate"), &UnderGenWorld3D::generate);
    ClassDB::bind_method(D_METHOD("cancel_generation"), &UnderGenWorld3D::cancel_generation);
    ClassDB::bind_method(D_METHOD("get_is_generating"), &UnderGenWorld3D::get_is_generating);

    // Callbacks bound for call_deferred safety
    ClassDB::bind_method(D_METHOD("_on_layout_completed", "outputs"), &UnderGenWorld3D::_on_layout_completed);
    ClassDB::bind_method(D_METHOD("_on_meshing_completed", "outputs"), &UnderGenWorld3D::_on_meshing_completed);
    ClassDB::bind_method(D_METHOD("_on_spawning_completed"), &UnderGenWorld3D::_on_spawning_completed);
    ClassDB::bind_method(D_METHOD("_on_generation_failed", "reason"), &UnderGenWorld3D::_on_generation_failed);

    // Inspector Properties
    ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "pipeline", PROPERTY_HINT_RESOURCE_TYPE, "UnderGenPipeline"), "set_pipeline", "get_pipeline");
    ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "grammar", PROPERTY_HINT_RESOURCE_TYPE, "Resource"), "set_grammar_override", "get_grammar_override");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "seed"), "set_generation_seed", "get_generation_seed");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "voxel_size"), "set_voxel_size", "get_voxel_size");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "surface_threshold"), "set_surface_threshold", "get_surface_threshold");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "generate_on_ready"), "set_generate_on_ready", "get_generate_on_ready");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "spawn_on_generation_complete"), "set_spawn_on_generation_complete", "get_spawn_on_generation_complete");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "cast_shadows"), "set_cast_shadows", "get_cast_shadows");
    ADD_PROPERTY(PropertyInfo(Variant::NODE_PATH, "parent_node"), "set_parent_node", "get_parent_node");
    ADD_PROPERTY(PropertyInfo(Variant::NODE_PATH, "multiplayer_spawner"), "set_multiplayer_spawner", "get_multiplayer_spawner");

    // Inspector button — toggle on to trigger generation
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "Generate Level"), "set_trigger_generate", "get_trigger_generate");

    // Debug Visualization
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "debug_show_zone_labels"), "set_debug_show_zone_labels", "get_debug_show_zone_labels");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "debug_zone_label_font_size", PROPERTY_HINT_RANGE, "8,128,1"), "set_debug_zone_label_font_size", "get_debug_zone_label_font_size");
    ADD_PROPERTY(PropertyInfo(Variant::COLOR, "debug_zone_label_color"), "set_debug_zone_label_color", "get_debug_zone_label_color");

    ClassDB::bind_method(D_METHOD("get_density_grid"), &UnderGenWorld3D::get_density_grid);
    ClassDB::bind_method(D_METHOD("get_semantic_graph"), &UnderGenWorld3D::get_semantic_graph);
    ClassDB::bind_method(D_METHOD("get_embedded_layout"), &UnderGenWorld3D::get_embedded_layout);
    ClassDB::bind_method(D_METHOD("get_geometry_plan"), &UnderGenWorld3D::get_geometry_plan);
    ClassDB::bind_method(D_METHOD("move_embedded_space", "id", "position", "rebuild"), &UnderGenWorld3D::move_embedded_space, DEFVAL(true));
    ClassDB::bind_method(D_METHOD("set_embedded_space_elevation", "id", "elevation", "rebuild"), &UnderGenWorld3D::set_embedded_space_elevation, DEFVAL(true));
    ClassDB::bind_method(D_METHOD("move_elevation_band", "band", "delta_y", "rebuild", "include_structural_spaces"), &UnderGenWorld3D::move_elevation_band, DEFVAL(true), DEFVAL(false));
    ClassDB::bind_method(D_METHOD("rebuild_spatial_geometry"), &UnderGenWorld3D::rebuild_spatial_geometry);
    ClassDB::bind_method(D_METHOD("get_last_context"), &UnderGenWorld3D::get_last_context);

    // Signals
    ADD_SIGNAL(MethodInfo("generation_started"));
    ADD_SIGNAL(MethodInfo("layout_completed"));
    ADD_SIGNAL(MethodInfo("meshing_completed"));
    ADD_SIGNAL(MethodInfo("spawning_completed"));
    ADD_SIGNAL(MethodInfo("generation_failed", PropertyInfo(Variant::STRING, "reason")));
    ADD_SIGNAL(MethodInfo("spatial_layout_changed", PropertyInfo(Variant::STRING, "source_id"), PropertyInfo(Variant::ARRAY, "dirty_regions"), PropertyInfo(Variant::INT, "revision")));
    ADD_SIGNAL(MethodInfo("spatial_geometry_rebuilt", PropertyInfo(Variant::ARRAY, "dirty_regions"), PropertyInfo(Variant::INT, "revision")));
}

void UnderGenWorld3D::set_pipeline(const Ref<UnderGenPipeline> &p_pipeline) {
    pipeline = p_pipeline;
}

Ref<UnderGenPipeline> UnderGenWorld3D::get_pipeline() const {
    return pipeline;
}

void UnderGenWorld3D::set_grammar_override(const Ref<Resource> &p_grammar) {
    grammar_override = p_grammar;
}

Ref<Resource> UnderGenWorld3D::get_grammar_override() const {
    return grammar_override;
}

void UnderGenWorld3D::set_generation_seed(int64_t p_seed) {
    generation_seed = p_seed;
}

int64_t UnderGenWorld3D::get_generation_seed() const {
    return generation_seed;
}

void UnderGenWorld3D::set_voxel_size(float p_size) {
    voxel_size = p_size;
}

float UnderGenWorld3D::get_voxel_size() const {
    return voxel_size;
}

void UnderGenWorld3D::set_surface_threshold(float p_threshold) {
    surface_threshold = p_threshold;
}

float UnderGenWorld3D::get_surface_threshold() const {
    return surface_threshold;
}

void UnderGenWorld3D::set_generate_on_ready(bool p_enabled) {
    generate_on_ready = p_enabled;
}

bool UnderGenWorld3D::get_generate_on_ready() const {
    return generate_on_ready;
}

void UnderGenWorld3D::set_spawn_on_generation_complete(bool p_enabled) {
    spawn_on_generation_complete = p_enabled;
}

bool UnderGenWorld3D::get_spawn_on_generation_complete() const {
    return spawn_on_generation_complete;
}

void UnderGenWorld3D::set_cast_shadows(bool p_cast_shadows) {
    cast_shadows = p_cast_shadows;
}

bool UnderGenWorld3D::get_cast_shadows() const {
    return cast_shadows;
}

void UnderGenWorld3D::set_parent_node(const NodePath &p_path) {
    parent_node_path = p_path;
}

NodePath UnderGenWorld3D::get_parent_node() const {
    return parent_node_path;
}

void UnderGenWorld3D::set_multiplayer_spawner(const NodePath &p_path) {
    multiplayer_spawner_path = p_path;
}

NodePath UnderGenWorld3D::get_multiplayer_spawner() const {
    return multiplayer_spawner_path;
}

// ── Inspector "button" ─────────────────────────────────────────────────

void UnderGenWorld3D::set_trigger_generate(bool p_val) {
    if (p_val) {
        generate();
    }
    _trigger_generate = false; // always reset — acts as a button, not a state
}

bool UnderGenWorld3D::get_trigger_generate() const {
    return false; // never shows as checked
}

// ── Debug Visualization ────────────────────────────────────────────────

void UnderGenWorld3D::set_debug_show_zone_labels(bool p_enabled) {
    debug_show_zone_labels = p_enabled;
}

bool UnderGenWorld3D::get_debug_show_zone_labels() const {
    return debug_show_zone_labels;
}

void UnderGenWorld3D::set_debug_zone_label_font_size(int p_size) {
    debug_zone_label_font_size = (p_size < 8) ? 8 : ((p_size > 128) ? 128 : p_size);
}

int UnderGenWorld3D::get_debug_zone_label_font_size() const {
    return debug_zone_label_font_size;
}

void UnderGenWorld3D::set_debug_zone_label_color(const Color &p_color) {
    debug_zone_label_color = p_color;
}

Color UnderGenWorld3D::get_debug_zone_label_color() const {
    return debug_zone_label_color;
}

void UnderGenWorld3D::generate() {
    if (pipeline.is_null()) {
        UtilityFunctions::printerr("UnderGenWorld3D: No pipeline assigned.");
        emit_signal("generation_failed", "No pipeline assigned.");
        return;
    }

    std::lock_guard<std::mutex> lock(thread_mutex);
    if (is_generating) {
        UtilityFunctions::print("UnderGenWorld3D: Generation already in progress. Cancelling and restarting.");
        cancel_generation();
    }

    // 1. Clear previous children
    for (int i = get_child_count() - 1; i >= 0; --i) {
        Node *child = get_child(i);
        if (child) {
            child->queue_free();
        }
    }

    is_generating = true;
    emit_signal("generation_started");

    // 2. Start generation thread
    gen_thread = std::thread(&UnderGenWorld3D::_run_generation_async, this, generation_seed);
}

void UnderGenWorld3D::cancel_generation() {
    // Note: Since C++ std::thread doesn't have cooperative cancel,
    // we let it complete or detach it. For safety in GDExtension,
    // we wait for it to join if it's running.
    if (gen_thread.joinable()) {
        gen_thread.join();
    }
    is_generating = false;
}

void UnderGenWorld3D::_run_generation_async(int64_t p_seed) {
    // ── Grammar override ──────────────────────────────────────────────────
    // If the user assigned a grammar resource on this node, inject it into
    // the pipeline's first Grammar Expander node so it supersedes whatever
    // .tres path was baked into the pipeline JSON.
    if (grammar_override.is_valid()) {
        Array nodes = pipeline->get_nodes();
        for (int i = 0; i < nodes.size(); ++i) {
            Ref<UnderGenNode> node = nodes[i];
            if (node.is_null()) continue;
            UnderGenGrammarNode *grammar_node =
                Object::cast_to<UnderGenGrammarNode>(node.ptr());
            if (grammar_node) {
                String path = grammar_override->get_path();
                if (!path.is_empty()) {
                    grammar_node->set_grammar_resource_path(path);
                    UtilityFunctions::print(
                        "UnderGenWorld3D: Grammar override → ", path);
                }
                break; // only override the first grammar node
            }
        }
    }

    // Propagate settings (surface threshold and voxel size) to all nodes in the pipeline
    if (pipeline.is_valid()) {
        Array nodes = pipeline->get_nodes();
        for (int i = 0; i < nodes.size(); ++i) {
            Ref<UnderGenNode> node = nodes[i];
            if (node.is_null()) continue;
            node->set("surface_threshold", surface_threshold);
            node->set("voxel_size", voxel_size);
            node->set("cast_shadows", cast_shadows);
        }
    }

    // Heavy processing runs on this background thread.
    // Inputs dictionary can contain global settings like seed and voxel_size.
    Dictionary initial_inputs;
    initial_inputs[0] = p_seed; // Seed mapped to Port 0 by convention
    initial_inputs[1] = voxel_size; // Voxel size mapped to Port 1 by convention

    Dictionary outputs;
    bool success = pipeline->execute_pipeline(initial_inputs, outputs);

    if (success) {
        // Deferred layout completed callback (transition to main thread)
        call_deferred("_on_layout_completed", outputs);
    } else {
        call_deferred("_on_generation_failed", "Pipeline execution failed.");
    }
}

void UnderGenWorld3D::_on_layout_completed(const Dictionary &outputs) {
    emit_signal("layout_completed");
    
    // Execute all Mesher nodes on the main thread
    if (pipeline.is_valid()) {
        Array nodes = pipeline->get_nodes();
        for (int i = 0; i < nodes.size(); ++i) {
            Ref<UnderGenNode> node = nodes[i];
            if (node.is_null()) continue;

            UnderGenMesherNode* mesher = Object::cast_to<UnderGenMesherNode>(node.ptr());
            if (mesher) {
                mesher->set_cast_shadows(cast_shadows);
                Dictionary inputs = pipeline->get_node_inputs(mesher->get_name());
                
                // Capture context and the mesher's own voxel_size for debug labels.
                // (The world node's voxel_size may differ from the mesher's,
                //  and the mesher places chunks using its own voxel_size.)
                Dictionary context = inputs.get(0, Dictionary());
                if (!context.is_empty()) {
                    _last_context = context;
                    _last_voxel_size = mesher->get_voxel_size();
                    update_gizmos();
                }

                Dictionary node_outputs;
                mesher->execute_with_parent(inputs, node_outputs, this);
            }
        }
    }
    
    call_deferred("_on_meshing_completed", outputs);
}

void UnderGenWorld3D::_on_meshing_completed(const Dictionary &outputs) {
    emit_signal("meshing_completed");

    Ref<DensityGrid> grid;
    Dictionary material_thicknesses;
    if (!_last_context.is_empty()) {
        grid = _last_context.get("grid", Ref<DensityGrid>());
        material_thicknesses = _last_context.get("material_thicknesses", Dictionary());
    }

    // Spawn debug zone labels if enabled
    if (debug_show_zone_labels) {
        _spawn_debug_zone_labels(_last_context);
    }

    // Prefer the final context consumed by the mesher. This lets downstream
    // semantic marker nodes replace authored route markers without older
    // intermediate node outputs reintroducing duplicates.
    vox_spawns.clear();
    if (_last_context.has("vox_spawns")) {
        vox_spawns = ((Array)_last_context["vox_spawns"]).duplicate();
    } else if (pipeline.is_valid()) {
        Array nodes = pipeline->get_nodes();
        for (int i = 0; i < nodes.size(); ++i) {
            Ref<UnderGenNode> node = nodes[i];
            if (node.is_null()) continue;
            Dictionary node_outputs = pipeline->get_node_outputs(node->get_name());
            Dictionary context = node_outputs.get(0, Dictionary());
            if (context.has("vox_spawns")) {
                Array spawns = context["vox_spawns"];
                for (int s_idx = 0; s_idx < spawns.size(); ++s_idx) {
                    Dictionary s_dict = spawns[s_idx];
                    if (!vox_spawns.has(s_dict)) {
                        vox_spawns.append(s_dict);
                    }
                }
            }
        }
    }

    // Adjust spawn positions to snap onto surfaces
    if (grid.is_valid()) {
        float surface = grid->get_surface_threshold();
        
        Vector3i directions[] = {
            Vector3i(0, -1, 0), // Down (floor)
            Vector3i(-1, 0, 0), // Left (wall)
            Vector3i(1, 0, 0),  // Right (wall)
            Vector3i(0, 0, -1), // Forward (wall)
            Vector3i(0, 0, 1),  // Backward (wall)
            Vector3i(0, 1, 0)   // Up (ceiling)
        };

        for (int i = 0; i < vox_spawns.size(); ++i) {
            Dictionary spawn_d = vox_spawns[i];
            Vector3 spawn_pos_grid = spawn_d.get("position", Vector3());
            Vector3i pos_spawn = Vector3i(Math::round(spawn_pos_grid.x), Math::round(spawn_pos_grid.y), Math::round(spawn_pos_grid.z));
            
            bool adjusted = false;
            for (const auto& direction : directions) {
                Vector3i pos_solid = pos_spawn + direction;
                if (grid->is_valid_position(pos_solid)) {
                    float val_solid = grid->get_cell(pos_solid);
                    if (val_solid > surface) {
                        float val_empty = grid->get_cell(pos_spawn);
                        float denominator = val_empty - val_solid;
                        float t = 0.5f;
                        if (Math::abs(denominator) > CMP_EPSILON) {
                            t = (surface - val_solid) / denominator;
                            t = Math::clamp(t, 0.0f, 1.0f);
                        }
                        
                        int mat = grid->get_material_id(pos_solid);
                        float thickness = 1.0f;
                        Variant key = mat;
                        if (material_thicknesses.has(key)) {
                            thickness = (float)material_thicknesses[key];
                        }
                        
                        float t_new = t * thickness;
                        
                        Vector3 pos_solid_f = Vector3(pos_solid);
                        Vector3 direction_f = Vector3(direction);
                        Vector3 adjusted_pos_grid = pos_solid_f - direction_f * t_new;
                        
                        spawn_d["position"] = adjusted_pos_grid * voxel_size;
                        adjusted = true;
                        break;
                    }
                }
            }
            if (!adjusted) {
                spawn_d["position"] = spawn_pos_grid * voxel_size;
            }
        }
    } else {
        for (int i = 0; i < vox_spawns.size(); ++i) {
            Dictionary spawn_d = vox_spawns[i];
            Vector3 spawn_pos_grid = spawn_d.get("position", Vector3());
            spawn_d["position"] = spawn_pos_grid * voxel_size;
        }
    }
    
    // Execute all Spawner nodes on the main thread if auto-spawn is enabled
    if (spawn_on_generation_complete) {
        spawn_scenes(this);
    }

    call_deferred("_on_spawning_completed");
}

void UnderGenWorld3D::_on_spawning_completed() {
    is_generating = false;
    
    // Join thread safely since it has finished executing
    if (gen_thread.joinable()) {
        gen_thread.detach(); // Detach since we are executing on the main thread callback
    }
    
    emit_signal("spawning_completed");
    UtilityFunctions::print("UnderGenWorld3D: Level generation completed successfully!");
}

void UnderGenWorld3D::_on_generation_failed(const String &reason) {
    is_generating = false;
    if (gen_thread.joinable()) {
        gen_thread.detach();
    }
    emit_signal("generation_failed", reason);
}

// ── Debug Zone Labels ──────────────────────────────────────────────────

void UnderGenWorld3D::_clear_debug_labels() {
    // Remove any existing debug label nodes (prefixed with "DebugZone_")
    for (int i = get_child_count() - 1; i >= 0; --i) {
        Node *child = get_child(i);
        if (child && child->get_name().begins_with("DebugZone_")) {
            child->queue_free();
        }
    }
}

void UnderGenWorld3D::_spawn_debug_zone_labels(const Dictionary &context) {
    if (context.is_empty()) {
        UtilityFunctions::printerr("UnderGenWorld3D: Cannot spawn zone labels — context is empty.");
        return;
    }

    _clear_debug_labels();

    Ref<DensityGrid> grid = context.get("grid", Ref<DensityGrid>());
    if (grid.is_null()) {
        UtilityFunctions::printerr("UnderGenWorld3D: Cannot spawn zone labels — no grid in context.");
        return;
    }

    Array rooms = context.get("rooms", Array());
    Array edges = context.get("edges", Array());
    float vs = _last_voxel_size;
    int label_idx = 0;

    // Helper lambda: create a Label3D with standard settings
    auto make_label = [&](const String &name_suffix, const String &text,
                          const Vector3 &world_pos, const Color &color) -> Label3D* {
        Label3D *lbl = memnew(Label3D);
        lbl->set_name("DebugZone_" + String::num_int64(label_idx++) + "_" + name_suffix);
        lbl->set_text(text);
        lbl->set_font_size(debug_zone_label_font_size);
        lbl->set_billboard_mode(BaseMaterial3D::BILLBOARD_ENABLED);
        lbl->set_outline_size(2);
        lbl->set_outline_modulate(Color(0.0f, 0.0f, 0.0f, 1.0f));
        lbl->set_modulate(color);
        lbl->set_position(world_pos);
        lbl->set_draw_flag(Label3D::FLAG_DISABLE_DEPTH_TEST, true);
        lbl->set_render_priority(127);
        return lbl;
    };

    // ── Room labels ─────────────────────────────────────────────────────
    // Build id→center lookup for edge midpoint calculation
    std::map<String, Vector3> room_centers;
    for (int i = 0; i < rooms.size(); ++i) {
        Dictionary room = rooms[i];
        if (room.is_empty()) continue;
        String rid = room.get("id", "");
        Vector3 center = room.get("center", Vector3());
        if (!rid.is_empty()) room_centers[rid] = center;
    }

    for (int i = 0; i < rooms.size(); ++i) {
        Dictionary room = rooms[i];
        if (room.is_empty()) continue;

        String room_type = room.get("type", "unknown");
        Vector3 center = room.get("center", Vector3());
        String room_id = room.get("id", "");

        String display_text = room_type;
        if (!room_id.is_empty()) display_text = room_type + "\n(" + room_id + ")";

        Vector3 world_center = center * vs + Vector3(0.0f, vs * 0.5f, 0.0f);
        Label3D *label = make_label(room_type, display_text, world_center, debug_zone_label_color);
        add_child(label);
    }

    // ── Edge (corridor) labels ──────────────────────────────────────────
    // Use a distinct cyan-green color for edges so they stand out from rooms
    Color edge_color = Color(0.2f, 1.0f, 0.8f, 1.0f);
    for (int i = 0; i < edges.size(); ++i) {
        Dictionary edge = edges[i];
        if (edge.is_empty()) continue;

        String edge_type = edge.get("type", "corridor");
        String from_id = edge.get("from", "");
        String to_id   = edge.get("to", "");

        // Find the midpoint between the two connected rooms
        auto fit = room_centers.find(from_id);
        auto tit = room_centers.find(to_id);
        if (fit == room_centers.end() || tit == room_centers.end()) continue;

        Vector3 midpoint = (fit->second + tit->second) * 0.5f;
        Vector3 world_pos = midpoint * vs + Vector3(0.0f, vs * 0.5f, 0.0f);

        String display = edge_type;
        if (from_id != to_id) display = edge_type + "\n" + from_id + " → " + to_id;

        Label3D *label = make_label("edge_" + edge_type, display, world_pos, edge_color);
        add_child(label);
    }

    int edge_count = (int)edges.size();
    int total = (int)rooms.size() + edge_count;
    UtilityFunctions::print("UnderGenWorld3D: Spawned ", rooms.size(), " room labels + ", edge_count, " edge labels (", total, " total).");
}

void UnderGenWorld3D::spawn_scenes(Node *parent_node) {
    if (!parent_node) {
        parent_node = this;
    }
    Node3D *parent_3d = Object::cast_to<Node3D>(parent_node);
    if (!parent_3d) {
        UtilityFunctions::printerr("UnderGenWorld3D::spawn_scenes: parent_node is not a Node3D.");
        return;
    }

    Node *target_parent = parent_node;
    if (!parent_node_path.is_empty()) {
        Node *custom_parent = parent_node->get_node_or_null(parent_node_path);
        if (custom_parent) {
            target_parent = custom_parent;
        } else {
            UtilityFunctions::printerr("UnderGenWorld3D: Custom parent node not found at path: ", parent_node_path);
        }
    }

    MultiplayerSpawner *spawner = nullptr;
    if (!multiplayer_spawner_path.is_empty()) {
        Node *spawner_node = parent_node->get_node_or_null(multiplayer_spawner_path);
        if (spawner_node) {
            spawner = Object::cast_to<MultiplayerSpawner>(spawner_node);
            if (!spawner) {
                UtilityFunctions::printerr("UnderGenWorld3D: Node at multiplayer_spawner path is not a MultiplayerSpawner.");
            }
        } else {
            UtilityFunctions::printerr("UnderGenWorld3D: MultiplayerSpawner not found at path: ", multiplayer_spawner_path);
        }
    }

    if (pipeline.is_valid()) {
        Array nodes = pipeline->get_nodes();
        for (int i = 0; i < nodes.size(); ++i) {
            Ref<UnderGenNode> node = nodes[i];
            if (node.is_null()) continue;

            UnderGenSceneSpawnerNode* spawner_node = Object::cast_to<UnderGenSceneSpawnerNode>(node.ptr());
            if (spawner_node) {
                Dictionary inputs = pipeline->get_node_inputs(spawner_node->get_name());
                Dictionary node_outputs;
                spawner_node->execute_with_parent(inputs, node_outputs, target_parent, spawner);
            }

            UnderGenMeshSpawnerNode* mesh_spawner_node = Object::cast_to<UnderGenMeshSpawnerNode>(node.ptr());
            if (mesh_spawner_node) {
                Dictionary inputs = pipeline->get_node_inputs(mesh_spawner_node->get_name());
                Dictionary node_outputs;
                mesh_spawner_node->execute_with_parent(inputs, node_outputs, target_parent);
            }
        }
    }
}

void UnderGenWorld3D::spawn_scenes_for_node(const String &node_name, Node *parent_node) {
    if (!parent_node) {
        parent_node = this;
    }
    Node3D *parent_3d = Object::cast_to<Node3D>(parent_node);
    if (!parent_3d) {
        UtilityFunctions::printerr("UnderGenWorld3D::spawn_scenes_for_node: parent_node is not a Node3D.");
        return;
    }

    Node *target_parent = parent_node;
    if (!parent_node_path.is_empty()) {
        Node *custom_parent = parent_node->get_node_or_null(parent_node_path);
        if (custom_parent) {
            target_parent = custom_parent;
        } else {
            UtilityFunctions::printerr("UnderGenWorld3D: Custom parent node not found at path: ", parent_node_path);
        }
    }

    MultiplayerSpawner *spawner = nullptr;
    if (!multiplayer_spawner_path.is_empty()) {
        Node *spawner_node = parent_node->get_node_or_null(multiplayer_spawner_path);
        if (spawner_node) {
            spawner = Object::cast_to<MultiplayerSpawner>(spawner_node);
            if (!spawner) {
                UtilityFunctions::printerr("UnderGenWorld3D: Node at multiplayer_spawner path is not a MultiplayerSpawner.");
            }
        } else {
            UtilityFunctions::printerr("UnderGenWorld3D: MultiplayerSpawner not found at path: ", multiplayer_spawner_path);
        }
    }

    if (pipeline.is_valid()) {
        Array nodes = pipeline->get_nodes();
        for (int i = 0; i < nodes.size(); ++i) {
            Ref<UnderGenNode> node = nodes[i];
            if (node.is_null()) continue;
            if (node->get_name() == node_name) {
                UnderGenSceneSpawnerNode* spawner_node = Object::cast_to<UnderGenSceneSpawnerNode>(node.ptr());
                if (spawner_node) {
                    Dictionary inputs = pipeline->get_node_inputs(node_name);
                    Dictionary node_outputs;
                    spawner_node->execute_with_parent(inputs, node_outputs, target_parent, spawner);
                } else {
                    UnderGenMeshSpawnerNode* mesh_spawner_node = Object::cast_to<UnderGenMeshSpawnerNode>(node.ptr());
                    if (mesh_spawner_node) {
                        Dictionary inputs = pipeline->get_node_inputs(node_name);
                        Dictionary node_outputs;
                        mesh_spawner_node->execute_with_parent(inputs, node_outputs, target_parent);
                    } else {
                        UtilityFunctions::printerr("UnderGenWorld3D::spawn_scenes_for_node: Node '", node_name, "' is not a supported spawner node.");
                    }
                }
                return;
            }
        }
        UtilityFunctions::printerr("UnderGenWorld3D::spawn_scenes_for_node: Node '", node_name, "' not found in pipeline.");
    }
}

Ref<UnderGenPointSet> UnderGenWorld3D::get_point_set_from_node(const String &node_name) const {
    if (pipeline.is_null()) {
        return Ref<UnderGenPointSet>();
    }
    Dictionary outputs = pipeline->get_node_outputs(node_name);
    Array keys = outputs.keys();
    for (int i = 0; i < keys.size(); i++) {
        Variant val = outputs[keys[i]];
        if (val.get_type() == Variant::OBJECT) {
            Ref<UnderGenPointSet> ps = val;
            if (ps.is_valid()) {
                return ps;
            }
        }
    }
    return Ref<UnderGenPointSet>();
}

Array UnderGenWorld3D::get_vox_spawns() const {
    return vox_spawns;
}

Ref<DensityGrid> UnderGenWorld3D::get_density_grid() const {
    if (_last_context.has("grid")) {
        return _last_context["grid"];
    }
    return Ref<DensityGrid>();
}

Ref<UnderGenSemanticGraph> UnderGenWorld3D::get_semantic_graph() const {
    return _last_context.get("semantic_graph", Ref<UnderGenSemanticGraph>());
}

Ref<UnderGenEmbeddedLayout> UnderGenWorld3D::get_embedded_layout() const {
    return _last_context.get("embedded_layout", Ref<UnderGenEmbeddedLayout>());
}

Ref<UnderGenGeometryPlan> UnderGenWorld3D::get_geometry_plan() const {
    return _last_context.get("geometry_plan", Ref<UnderGenGeometryPlan>());
}

bool UnderGenWorld3D::move_embedded_space(const String &p_id, const Vector3 &p_position, bool p_rebuild) {
    if (is_generating) return false;
    Ref<UnderGenEmbeddedLayout> layout = get_embedded_layout();
    if (layout.is_null() || !layout->move_space(p_id, p_position)) return false;
    layout->validate_layout();
    emit_signal("spatial_layout_changed", p_id, layout->get_dirty_regions(), layout->get_revision());
    return !p_rebuild || rebuild_spatial_geometry();
}

bool UnderGenWorld3D::set_embedded_space_elevation(const String &p_id, float p_elevation, bool p_rebuild) {
    Ref<UnderGenEmbeddedLayout> layout = get_embedded_layout();
    if (layout.is_null()) return false;
    Ref<UnderGenEmbeddedSpace> space = layout->find_space(p_id);
    if (space.is_null()) return false;
    Vector3 position = space->get_position(); position.y = p_elevation;
    return move_embedded_space(p_id, position, p_rebuild);
}

bool UnderGenWorld3D::move_elevation_band(int p_band, float p_delta_y, bool p_rebuild, bool p_include_structural_spaces) {
    if (is_generating) return false;
    Ref<UnderGenEmbeddedLayout> layout = get_embedded_layout();
    if (layout.is_null() || !layout->move_elevation_band(p_band, p_delta_y, p_include_structural_spaces)) return false;
    layout->validate_layout();
    emit_signal("spatial_layout_changed", "band:" + String::num_int64(p_band), layout->get_dirty_regions(), layout->get_revision());
    return !p_rebuild || rebuild_spatial_geometry();
}

bool UnderGenWorld3D::rebuild_spatial_geometry() {
    if (is_generating || pipeline.is_null()) return false;
    Ref<UnderGenEmbeddedLayout> layout = get_embedded_layout();
    if (layout.is_null()) return false;
    UnderGenGeometryPlannerNode *planner = nullptr;
    UnderGenGeometryRealizerNode *realizer = nullptr;
    UnderGenGameplayMarkerNode *gameplay_markers = nullptr;
    Array nodes = pipeline->get_nodes();
    for (int i = 0; i < nodes.size(); ++i) {
        Ref<UnderGenNode> node = nodes[i]; if (node.is_null()) continue;
        if (!planner) planner = Object::cast_to<UnderGenGeometryPlannerNode>(node.ptr());
        if (!realizer) realizer = Object::cast_to<UnderGenGeometryRealizerNode>(node.ptr());
        if (!gameplay_markers) gameplay_markers = Object::cast_to<UnderGenGameplayMarkerNode>(node.ptr());
    }
    if (!planner || !realizer) return false;
    Array dirty_regions = layout->get_dirty_regions();
    if (dirty_regions.is_empty()) dirty_regions.append(layout->get_bounds());
    Ref<UnderGenGeometryPlan> plan = planner->build_plan(layout);
    if (plan.is_null()) return false;
    _last_context = realizer->rebuild_dirty_regions(_last_context, plan, dirty_regions);
    if (gameplay_markers) {
        Dictionary marker_inputs;
        marker_inputs[0] = _last_context;
        Dictionary marker_outputs;
        gameplay_markers->execute(marker_inputs, marker_outputs);
        if (marker_outputs.has(0) && marker_outputs[0].get_type() == Variant::DICTIONARY) {
            _last_context = marker_outputs[0];
            vox_spawns = ((Array)_last_context.get("vox_spawns", Array())).duplicate();
        }
    }
    _remesh_dirty_chunks(dirty_regions);
    layout->clear_dirty_regions();
    update_gizmos();
    emit_signal("spatial_geometry_rebuilt", dirty_regions, layout->get_revision());
    return true;
}

void UnderGenWorld3D::_remesh_dirty_chunks(const Array &p_dirty_regions) {
    for (int i = 0; i < get_child_count(); ++i) {
        Node *child = get_child(i);
        Vector3i offset;
        int chunk_size = 0;
        MCChunk *mc = Object::cast_to<MCChunk>(child);
        DCChunk *dc = Object::cast_to<DCChunk>(child);
        if (mc) { offset = mc->get_chunk_grid_offset(); chunk_size = mc->get_chunk_size(); }
        else if (dc) { offset = dc->get_chunk_grid_offset(); chunk_size = dc->get_chunk_size(); }
        else continue;
        AABB chunk_bounds(Vector3(offset) - Vector3(1, 1, 1), Vector3(chunk_size + 2, chunk_size + 2, chunk_size + 2));
        bool touched = false;
        for (int r = 0; r < p_dirty_regions.size(); ++r) {
            if (p_dirty_regions[r].get_type() == Variant::AABB && chunk_bounds.intersects((AABB)p_dirty_regions[r])) { touched = true; break; }
        }
        if (touched) child->call_deferred("generate_mesh_from_density_grid");
    }
}

} // namespace godot
