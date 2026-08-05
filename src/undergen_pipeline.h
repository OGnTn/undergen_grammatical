#ifndef UNDERGEN_PIPELINE_H
#define UNDERGEN_PIPELINE_H

#include <godot_cpp/classes/resource.hpp>
#include <godot_cpp/variant/array.hpp>
#include <godot_cpp/variant/dictionary.hpp>
#include <godot_cpp/variant/string.hpp>
#include <godot_cpp/variant/typed_array.hpp>
#include <godot_cpp/variant/vector2.hpp>
#include "undergen_node.h"
#include <map>
#include <vector>

namespace godot {

class UnderGenPipeline : public Resource {
    GDCLASS(UnderGenPipeline, Resource);

private:
    Array nodes; // List of UnderGenNode resources
    Array connections; // Connection dictionaries: { "from_node": string, "from_port": int, "to_node": string, "to_port": int }

protected:
    static void _bind_methods();

public:
    UnderGenPipeline();
    virtual ~UnderGenPipeline();

    // ── Node list management ──────────────────────────────────────────────
    void set_nodes(const Array &p_nodes);
    Array get_nodes() const;
    void add_pipeline_node(const Ref<UnderGenNode> &p_node);
    void remove_pipeline_node(const Ref<UnderGenNode> &p_node);

    // ── Connections management ────────────────────────────────────────────
    void set_connections(const Array &p_connections);
    Array get_connections() const;
    void add_connection(const String &from_node, int from_port, const String &to_node, int to_port);
    void remove_connection(const String &from_node, int from_port, const String &to_node, int to_port);

    // ── Pipeline Execution Engine ─────────────────────────────────────────
    bool execute_pipeline(const Dictionary &initial_inputs, Dictionary &out_final_outputs);
    Dictionary execute(const Dictionary &initial_inputs);

    Dictionary get_node_inputs(const String &node_name) const;
    Dictionary get_node_outputs(const String &node_name) const;

    // ── Builder API ───────────────────────────────────────────────────────
    /// Create a node of a given class name, configure it via a dictionary of
    /// property→value entries, and add it to the pipeline.  Returns the node
    /// for further chaining.
    ///
    /// Usage (GDScript):
    ///   pipeline.add_node_of_type("UnderGenBSPPlacerNode", "BSP_Placer",
    ///       {"grid_size_x": 64, "grid_size_y": 20, "grid_size_z": 64},
    ///       Vector2(300, 200))
    Ref<UnderGenNode> add_node_of_type(const String &class_name,
                                        const String &node_name,
                                        const Dictionary &properties = Dictionary(),
                                        const Vector2 &editor_pos = Vector2());

    /// Quick connect two nodes by name.
    void connect(const String &from_name, int from_port,
                 const String &to_name, int to_port);

    // ── Serialization ─────────────────────────────────────────────────────
    /// Serialize the entire pipeline to a Godot Dictionary.
    Dictionary to_dictionary() const;

    /// Populate from a Godot Dictionary.
    /// WARNING: overwrites existing nodes/connections.
    void from_dictionary(const Dictionary &d);

    /// Serialize to a pretty-printed JSON string.
    String to_json(bool indent = true) const;

    /// Populate from a JSON string. Returns OK on success.
    Error from_json(const String &json_str,
                    int *r_err_line = nullptr,
                    String *r_err_message = nullptr);

    /// ClassDB-bindable overload (no pointer params).
    Error _from_json_bound(const String &json_str);

    /// Save this resource to disk as a .tres file.
    Error save_to_file(const String &path);

    /// Static factory: load pipeline from a .json file on disk.
    static Ref<UnderGenPipeline> load_from_json_file(const String &path);

private:
    Dictionary node_inputs_cache;
    Dictionary node_outputs_cache;

    // Helper to topologically sort the nodes
    bool _get_topological_order(std::vector<String> &r_order, const std::map<String, Ref<UnderGenNode>> &nodes_map);
};

} // namespace godot

#endif // UNDERGEN_PIPELINE_H
