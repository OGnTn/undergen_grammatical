#include "undergen_pipeline.h"
#include <godot_cpp/variant/utility_functions.hpp>
#include <queue>
#include <set>

namespace godot {

UnderGenPipeline::UnderGenPipeline() {}
UnderGenPipeline::~UnderGenPipeline() {}

void UnderGenPipeline::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_nodes", "nodes"), &UnderGenPipeline::set_nodes);
    ClassDB::bind_method(D_METHOD("get_nodes"), &UnderGenPipeline::get_nodes);
    ClassDB::bind_method(D_METHOD("add_pipeline_node", "node"), &UnderGenPipeline::add_pipeline_node);
    ClassDB::bind_method(D_METHOD("remove_pipeline_node", "node"), &UnderGenPipeline::remove_pipeline_node);

    ClassDB::bind_method(D_METHOD("set_connections", "connections"), &UnderGenPipeline::set_connections);
    ClassDB::bind_method(D_METHOD("get_connections"), &UnderGenPipeline::get_connections);
    ClassDB::bind_method(D_METHOD("add_connection", "from_node", "from_port", "to_node", "to_port"), &UnderGenPipeline::add_connection);
    ClassDB::bind_method(D_METHOD("remove_connection", "from_node", "from_port", "to_node", "to_port"), &UnderGenPipeline::remove_connection);

    ClassDB::bind_method(D_METHOD("get_node_inputs", "node_name"), &UnderGenPipeline::get_node_inputs);
    ClassDB::bind_method(D_METHOD("get_node_outputs", "node_name"), &UnderGenPipeline::get_node_outputs);

    ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "nodes"), "set_nodes", "get_nodes");
    ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "connections"), "set_connections", "get_connections");
}

void UnderGenPipeline::set_nodes(const Array &p_nodes) {
    nodes = p_nodes;
}

Array UnderGenPipeline::get_nodes() const {
    return nodes;
}

void UnderGenPipeline::add_pipeline_node(const Ref<UnderGenNode> &p_node) {
    if (p_node.is_valid() && !nodes.has(p_node)) {
        nodes.append(p_node);
    }
}

void UnderGenPipeline::remove_pipeline_node(const Ref<UnderGenNode> &p_node) {
    if (p_node.is_valid()) {
        nodes.erase(p_node);
        // Also remove any connections involving this node
        String node_name = p_node->get_name();
        Array filtered_conn;
        for (int i = 0; i < connections.size(); ++i) {
            Dictionary conn = connections[i];
            String from = conn.get("from_node", "");
            String to = conn.get("to_node", "");
            if (from != node_name && to != node_name) {
                filtered_conn.append(conn);
            }
        }
        connections = filtered_conn;
    }
}

void UnderGenPipeline::set_connections(const Array &p_connections) {
    connections = p_connections;
}

Array UnderGenPipeline::get_connections() const {
    return connections;
}

void UnderGenPipeline::add_connection(const String &from_node, int from_port, const String &to_node, int to_port) {
    // Check if it already exists
    for (int i = 0; i < connections.size(); ++i) {
        Dictionary conn = connections[i];
        if ((String)conn.get("from_node", "") == from_node &&
            (int)conn.get("from_port", 0) == from_port &&
            (String)conn.get("to_node", "") == to_node &&
            (int)conn.get("to_port", 0) == to_port) {
            return; // Already exists
        }
    }

    Dictionary new_conn;
    new_conn["from_node"] = from_node;
    new_conn["from_port"] = from_port;
    new_conn["to_node"] = to_node;
    new_conn["to_port"] = to_port;
    connections.append(new_conn);
}

void UnderGenPipeline::remove_connection(const String &from_node, int from_port, const String &to_node, int to_port) {
    for (int i = 0; i < connections.size(); ++i) {
        Dictionary conn = connections[i];
        if ((String)conn.get("from_node", "") == from_node &&
            (int)conn.get("from_port", 0) == from_port &&
            (String)conn.get("to_node", "") == to_node &&
            (int)conn.get("to_port", 0) == to_port) {
            connections.remove_at(i);
            break;
        }
    }
}

bool UnderGenPipeline::execute_pipeline(const Dictionary &initial_inputs, Dictionary &out_final_outputs) {
    // Clear caches
    node_inputs_cache.clear();
    node_outputs_cache.clear();

    // 1. Map names to Ref<UnderGenNode>
    std::map<String, Ref<UnderGenNode>> nodes_map;
    for (int i = 0; i < nodes.size(); ++i) {
        Ref<UnderGenNode> n = nodes[i];
        if (n.is_valid()) {
            nodes_map[n->get_name()] = n;
        }
    }

    // 2. Perform topological sort
    std::vector<String> execution_order;
    if (!_get_topological_order(execution_order, nodes_map)) {
        UtilityFunctions::printerr("UnderGenPipeline: Cycle detected or invalid connections in pipeline graph!");
        return false;
    }

    // 3. Setup outputs cache
    std::map<String, Dictionary> node_outputs;
    
    // Seed initial inputs to a virtual "input" node or directly to initial nodes
    node_outputs["input"] = initial_inputs.duplicate();

    // 4. Execute nodes sequentially
    for (const String &node_name : execution_order) {
        Ref<UnderGenNode> node = nodes_map[node_name];
        if (node.is_null()) {
            continue;
        }

        // Gather inputs for this node
        Dictionary inputs;
        for (int i = 0; i < connections.size(); ++i) {
            Dictionary conn = connections[i];
            String to = conn.get("to_node", "");
            if (to == node_name) {
                String from = conn.get("from_node", "");
                int from_port = conn.get("from_port", 0);
                int to_port = conn.get("to_port", 0);

                if (node_outputs.count(from)) {
                    inputs[to_port] = node_outputs[from].get(from_port, Variant());
                }
            }
        }

        // Cache the inputs for this node
        node_inputs_cache[node_name] = inputs.duplicate();

        // Execute the node
        Dictionary outputs;
        node->execute(inputs, outputs);

        // Cache and store outputs
        node_outputs_cache[node_name] = outputs.duplicate();
        node_outputs[node_name] = outputs;
    }

    // 5. Gather final outputs from endpoints (or nodes connected to nothing/output node)
    // If there is an output node, we look at its inputs, or we return the outputs of the last node
    if (!execution_order.empty()) {
        String last_node = execution_order.back();
        out_final_outputs = node_outputs[last_node];
    }

    return true;
}

Dictionary UnderGenPipeline::get_node_inputs(const String &node_name) const {
    return node_inputs_cache.get(node_name, Dictionary());
}

Dictionary UnderGenPipeline::get_node_outputs(const String &node_name) const {
    return node_outputs_cache.get(node_name, Dictionary());
}


bool UnderGenPipeline::_get_topological_order(std::vector<String> &r_order, const std::map<String, Ref<UnderGenNode>> &nodes_map) {
    std::map<String, std::vector<String>> adj;
    std::map<String, int> in_degree;

    // Initialize degrees
    for (const auto &pair : nodes_map) {
        in_degree[pair.first] = 0;
    }

    // Build adjacency list
    for (int i = 0; i < connections.size(); ++i) {
        Dictionary conn = connections[i];
        String from = conn.get("from_node", "");
        String to = conn.get("to_node", "");

        if (nodes_map.count(from) && nodes_map.count(to)) {
            adj[from].push_back(to);
            in_degree[to]++;
        }
    }

    // Process nodes with in-degree 0
    std::queue<String> q;
    for (const auto &pair : nodes_map) {
        if (in_degree[pair.first] == 0) {
            q.push(pair.first);
        }
    }

    while (!q.empty()) {
        String u = q.front();
        q.pop();
        r_order.push_back(u);

        for (const String &v : adj[u]) {
            in_degree[v]--;
            if (in_degree[v] == 0) {
                q.push(v);
            }
        }
    }

    // If order size is less than nodes size, a cycle exists
    return r_order.size() == nodes_map.size();
}

} // namespace godot
