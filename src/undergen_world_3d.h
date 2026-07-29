#ifndef UNDERGEN_WORLD_3D_H
#define UNDERGEN_WORLD_3D_H

#include <godot_cpp/classes/node3d.hpp>
#include <godot_cpp/classes/ref.hpp>
#include <godot_cpp/variant/node_path.hpp>
#include "undergen_pipeline.h"
#include <thread>
#include <mutex>

namespace godot {

class DensityGrid;
class UnderGenPointSet;

class UnderGenWorld3D : public Node3D {
    GDCLASS(UnderGenWorld3D, Node3D);

private:
    Ref<UnderGenPipeline> pipeline;
    Ref<Resource> grammar_override;   // if set, supersedes the pipeline's Grammar Expander node
    int64_t generation_seed = 12345;
    float voxel_size = 1.0f;
    float surface_threshold = 0.0f;
    bool generate_on_ready = false;
    bool spawn_on_generation_complete = true;
    bool cast_shadows = true;
    Array vox_spawns;
    NodePath parent_node_path;
    NodePath multiplayer_spawner_path;

    // Inspector "button" — set to true triggers generation
    bool _trigger_generate = false;

    // Debug Visualization
    bool debug_show_zone_labels = false;
    int debug_zone_label_font_size = 24;
    Color debug_zone_label_color = Color(1.0f, 1.0f, 0.0f, 1.0f); // yellow
    Dictionary _last_context; // stored for debug label spawning
    float _last_voxel_size = 1.0f; // captured from mesher during layout

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

    void set_grammar_override(const Ref<Resource> &p_grammar);
    Ref<Resource> get_grammar_override() const;

    void set_generation_seed(int64_t p_seed);
    int64_t get_generation_seed() const;

    void set_voxel_size(float p_size);
    float get_voxel_size() const;

    void set_surface_threshold(float p_threshold);
    float get_surface_threshold() const;

    void set_generate_on_ready(bool p_enabled);
    bool get_generate_on_ready() const;

    void set_spawn_on_generation_complete(bool p_enabled);
    bool get_spawn_on_generation_complete() const;

    void set_cast_shadows(bool p_cast_shadows);
    bool get_cast_shadows() const;

    void set_parent_node(const NodePath &p_path);
    NodePath get_parent_node() const;

    void set_multiplayer_spawner(const NodePath &p_path);
    NodePath get_multiplayer_spawner() const;

    void spawn_scenes(Node *parent_node = nullptr);
    void spawn_scenes_for_node(const String &node_name, Node *parent_node = nullptr);
    Ref<UnderGenPointSet> get_point_set_from_node(const String &node_name) const;
    Array get_vox_spawns() const;

    // Inspector "button" — set true to trigger generation
    void set_trigger_generate(bool p_val);
    bool get_trigger_generate() const;

    // Debug Visualization
    void set_debug_show_zone_labels(bool p_enabled);
    bool get_debug_show_zone_labels() const;
    void set_debug_zone_label_font_size(int p_size);
    int get_debug_zone_label_font_size() const;
    void set_debug_zone_label_color(const Color &p_color);
    Color get_debug_zone_label_color() const;

    // Generation entry point
    void generate();
    void cancel_generation();
    bool get_is_generating() const { return is_generating; }

    Ref<DensityGrid> get_density_grid() const;
    Dictionary get_last_context() const { return _last_context; }

private:
    // Run on background thread
    void _run_generation_async(int64_t p_seed);

    // Main thread callbacks (called via call_deferred)
    void _on_layout_completed(const Dictionary &outputs);
    void _on_meshing_completed(const Dictionary &outputs);
    void _on_spawning_completed();
    void _on_generation_failed(const String &reason);

    // Debug helpers
    void _spawn_debug_zone_labels(const Dictionary &context);
    void _clear_debug_labels();
};

} // namespace godot

#endif // UNDERGEN_WORLD_3D_H
