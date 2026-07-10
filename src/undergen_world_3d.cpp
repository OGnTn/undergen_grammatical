#include "undergen_world_3d.h"
#include "undergen_mesher_node.h"
#include "undergen_scene_spawner_node.h"
#include "undergen_grammar_node.h"
#include "undergen_bsp_placer_node.h"
#include "density_grid.h"
#include "undergen_point_set.h"
#include "undergen_vox_stamp_node.h"
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
    ADD_PROPERTY(PropertyInfo(Variant::NODE_PATH, "parent_node"), "set_parent_node", "get_parent_node");
    ADD_PROPERTY(PropertyInfo(Variant::NODE_PATH, "multiplayer_spawner"), "set_multiplayer_spawner", "get_multiplayer_spawner");

    // Inspector button — toggle on to trigger generation
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "Generate Level"), "set_trigger_generate", "get_trigger_generate");

    // Debug Visualization
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "debug_show_zone_labels"), "set_debug_show_zone_labels", "get_debug_show_zone_labels");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "debug_zone_label_font_size", PROPERTY_HINT_RANGE, "8,128,1"), "set_debug_zone_label_font_size", "get_debug_zone_label_font_size");
    ADD_PROPERTY(PropertyInfo(Variant::COLOR, "debug_zone_label_color"), "set_debug_zone_label_color", "get_debug_zone_label_color");

    // Signals
    ADD_SIGNAL(MethodInfo("generation_started"));
    ADD_SIGNAL(MethodInfo("layout_completed"));
    ADD_SIGNAL(MethodInfo("meshing_completed"));
    ADD_SIGNAL(MethodInfo("spawning_completed"));
    ADD_SIGNAL(MethodInfo("generation_failed", PropertyInfo(Variant::STRING, "reason")));
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
                Dictionary inputs = pipeline->get_node_inputs(mesher->get_name());
                
                // Capture context and the mesher's own voxel_size for debug labels.
                // (The world node's voxel_size may differ from the mesher's,
                //  and the mesher places chunks using its own voxel_size.)
                Dictionary context = inputs.get(0, Dictionary());
                if (!context.is_empty()) {
                    _last_context = context;
                    _last_voxel_size = mesher->get_voxel_size();
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

    // Spawn debug zone labels if enabled
    if (debug_show_zone_labels) {
        _spawn_debug_zone_labels(_last_context);
        _last_context = Dictionary(); // Release reference safely without clearing shared dictionary data
    }

    // Gather all vox spawns from any VoxStamper nodes
    vox_spawns.clear();
    if (pipeline.is_valid()) {
        Array nodes = pipeline->get_nodes();
        for (int i = 0; i < nodes.size(); ++i) {
            Ref<UnderGenNode> node = nodes[i];
            if (node.is_null()) continue;
            UnderGenVoxStampNode* stamper = Object::cast_to<UnderGenVoxStampNode>(node.ptr());
            if (stamper) {
                Dictionary node_outputs = pipeline->get_node_outputs(stamper->get_name());
                Dictionary context = node_outputs.get(0, Dictionary());
                if (context.has("vox_spawns")) {
                    Array spawns = context["vox_spawns"];
                    vox_spawns.append_array(spawns);
                }
            }
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
                    UtilityFunctions::printerr("UnderGenWorld3D::spawn_scenes_for_node: Node '", node_name, "' is not an UnderGenSceneSpawnerNode.");
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

} // namespace godot
