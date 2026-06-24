#ifndef UNDERGEN_WORLD_3D_H
#define UNDERGEN_WORLD_3D_H

#include <godot_cpp/classes/node3d.hpp>
#include <godot_cpp/classes/ref.hpp>
#include "undergen_pipeline.h"
#include <thread>
#include <mutex>

namespace godot {

class UnderGenWorld3D : public Node3D {
    GDCLASS(UnderGenWorld3D, Node3D);

private:
    Ref<UnderGenPipeline> pipeline;
    int64_t generation_seed = 12345;
    float voxel_size = 1.0f;
    bool generate_on_ready = false;

    // Threading State
    std::thread gen_thread;
    std::mutex thread_mutex;
    bool is_generating = false;

protected:
    static void _bind_methods();
    void _notification(int p_what);

public:
    UnderGenWorld3D();
    virtual ~UnderGenWorld3D();

    // Getters / Setters
    void set_pipeline(const Ref<UnderGenPipeline> &p_pipeline);
    Ref<UnderGenPipeline> get_pipeline() const;

    void set_generation_seed(int64_t p_seed);
    int64_t get_generation_seed() const;

    void set_voxel_size(float p_size);
    float get_voxel_size() const;

    void set_generate_on_ready(bool p_enabled);
    bool get_generate_on_ready() const;

    // Generation entry point
    void generate();
    void cancel_generation();
    bool get_is_generating() const { return is_generating; }

private:
    // Run on background thread
    void _run_generation_async(int64_t p_seed);

    // Main thread callbacks (called via call_deferred)
    void _on_layout_completed(const Dictionary &outputs);
    void _on_meshing_completed(const Dictionary &outputs);
    void _on_spawning_completed();
    void _on_generation_failed(const String &reason);
};

} // namespace godot

#endif // UNDERGEN_WORLD_3D_H
