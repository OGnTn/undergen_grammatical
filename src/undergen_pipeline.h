#ifndef UNDERGEN_PIPELINE_H
#define UNDERGEN_PIPELINE_H

#include <godot_cpp/classes/resource.hpp>
#include <godot_cpp/variant/array.hpp>
#include <godot_cpp/variant/dictionary.hpp>
#include <godot_cpp/variant/string.hpp>
#include <godot_cpp/variant/typed_array.hpp>
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

    // Node list management
    void set_nodes(const Array &p_nodes);
    Array get_nodes() const;
    void add_pipeline_node(const Ref<UnderGenNode> &p_node);
    void remove_pipeline_node(const Ref<UnderGenNode> &p_node);

    // Connections management
    void set_connections(const Array &p_connections);
    Array get_connections() const;
    void add_connection(const String &from_node, int from_port, const String &to_node, int to_port);
    void remove_connection(const String &from_node, int from_port, const String &to_node, int to_port);

    // Pipeline Execution Engine
    bool execute_pipeline(const Dictionary &initial_inputs, Dictionary &out_final_outputs);

    Dictionary get_node_inputs(const String &node_name) const;
    Dictionary get_node_outputs(const String &node_name) const;

private:
    Dictionary node_inputs_cache;
    Dictionary node_outputs_cache;

    // Helper to topologically sort the nodes
    bool _get_topological_order(std::vector<String> &r_order, const std::map<String, Ref<UnderGenNode>> &nodes_map);
};

} // namespace godot

#endif // UNDERGEN_PIPELINE_H

