#include "undergen_world_3d.h"
#include "undergen_mesher_node.h"
#include "undergen_scene_spawner_node.h"
#include <godot_cpp/variant/utility_functions.hpp>

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

    ClassDB::bind_method(D_METHOD("set_generation_seed", "seed"), &UnderGenWorld3D::set_generation_seed);
    ClassDB::bind_method(D_METHOD("get_generation_seed"), &UnderGenWorld3D::get_generation_seed);

    ClassDB::bind_method(D_METHOD("set_voxel_size", "size"), &UnderGenWorld3D::set_voxel_size);
    ClassDB::bind_method(D_METHOD("get_voxel_size"), &UnderGenWorld3D::get_voxel_size);

    ClassDB::bind_method(D_METHOD("set_generate_on_ready", "enabled"), &UnderGenWorld3D::set_generate_on_ready);
    ClassDB::bind_method(D_METHOD("get_generate_on_ready"), &UnderGenWorld3D::get_generate_on_ready);

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
    ADD_PROPERTY(PropertyInfo(Variant::INT, "seed"), "set_generation_seed", "get_generation_seed");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "voxel_size"), "set_voxel_size", "get_voxel_size");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "generate_on_ready"), "set_generate_on_ready", "get_generate_on_ready");

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

void UnderGenWorld3D::set_generate_on_ready(bool p_enabled) {
    generate_on_ready = p_enabled;
}

bool UnderGenWorld3D::get_generate_on_ready() const {
    return generate_on_ready;
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
                Dictionary node_outputs;
                mesher->execute_with_parent(inputs, node_outputs, this);
            }
        }
    }
    
    call_deferred("_on_meshing_completed", outputs);
}

void UnderGenWorld3D::_on_meshing_completed(const Dictionary &outputs) {
    emit_signal("meshing_completed");
    
    // Execute all Spawner nodes on the main thread
    if (pipeline.is_valid()) {
        Array nodes = pipeline->get_nodes();
        for (int i = 0; i < nodes.size(); ++i) {
            Ref<UnderGenNode> node = nodes[i];
            if (node.is_null()) continue;

            UnderGenSceneSpawnerNode* spawner = Object::cast_to<UnderGenSceneSpawnerNode>(node.ptr());
            if (spawner) {
                Dictionary inputs = pipeline->get_node_inputs(spawner->get_name());
                Dictionary node_outputs;
                spawner->execute_with_parent(inputs, node_outputs, this);
            }
        }
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

} // namespace godot
