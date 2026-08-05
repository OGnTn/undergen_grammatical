#include "undergen_pipeline.h"
#include <godot_cpp/classes/json.hpp>
#include <godot_cpp/classes/file_access.hpp>
#include <godot_cpp/classes/resource_saver.hpp>
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
    ClassDB::bind_method(D_METHOD("execute", "initial_inputs"), &UnderGenPipeline::execute);

    // Builder
    ClassDB::bind_method(D_METHOD("add_node_of_type", "class_name", "node_name", "properties", "editor_pos"),
                         &UnderGenPipeline::add_node_of_type,
                         DEFVAL(Dictionary()), DEFVAL(Vector2()));
    ClassDB::bind_method(D_METHOD("connect", "from_name", "from_port", "to_name", "to_port"),
                         &UnderGenPipeline::connect);

    // Serialization
    ClassDB::bind_method(D_METHOD("to_dictionary"), &UnderGenPipeline::to_dictionary);
    ClassDB::bind_method(D_METHOD("from_dictionary", "d"), &UnderGenPipeline::from_dictionary);
    ClassDB::bind_method(D_METHOD("to_json", "indent"), &UnderGenPipeline::to_json, DEFVAL(true));
    ClassDB::bind_method(D_METHOD("from_json", "json_str"), &UnderGenPipeline::_from_json_bound);
    ClassDB::bind_method(D_METHOD("save_to_file", "path"), &UnderGenPipeline::save_to_file);
    ClassDB::bind_static_method("UnderGenPipeline", D_METHOD("load_from_json_file", "path"), &UnderGenPipeline::load_from_json_file);

    ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "nodes"), "set_nodes", "get_nodes");
    ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "connections"), "set_connections", "get_connections");
}

void UnderGenPipeline::set_nodes(const Array &p_nodes) {
    for (int i = 0; i < nodes.size(); ++i) {
        Ref<UnderGenNode> old_node = nodes[i];
        if (old_node.is_valid() && old_node->is_connected("changed", Callable(this, "emit_changed"))) {
            old_node->disconnect("changed", Callable(this, "emit_changed"));
        }
    }
    nodes = p_nodes;
    for (int i = 0; i < nodes.size(); ++i) {
        Ref<UnderGenNode> new_node = nodes[i];
        if (new_node.is_valid()) {
            new_node->connect("changed", Callable(this, "emit_changed"));
        }
    }
}

Array UnderGenPipeline::get_nodes() const {
    return nodes;
}

void UnderGenPipeline::add_pipeline_node(const Ref<UnderGenNode> &p_node) {
    if (p_node.is_valid() && !nodes.has(p_node)) {
        nodes.append(p_node);
        p_node->connect("changed", Callable(this, "emit_changed"));
    }
}

void UnderGenPipeline::remove_pipeline_node(const Ref<UnderGenNode> &p_node) {
    if (p_node.is_valid()) {
        nodes.erase(p_node);
        if (p_node->is_connected("changed", Callable(this, "emit_changed"))) {
            p_node->disconnect("changed", Callable(this, "emit_changed"));
        }
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
        Dictionary inputs = node->get_pipeline_input_defaults(initial_inputs);
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

Dictionary UnderGenPipeline::execute(const Dictionary &initial_inputs) {
    Dictionary outputs;
    if (!execute_pipeline(initial_inputs, outputs)) return Dictionary();
    return outputs;
}

Dictionary UnderGenPipeline::get_node_inputs(const String &node_name) const {
    return node_inputs_cache.get(node_name, Dictionary());
}

Dictionary UnderGenPipeline::get_node_outputs(const String &node_name) const {
    return node_outputs_cache.get(node_name, Dictionary());
}


// ── Builder API ───────────────────────────────────────────────────────────────

Ref<UnderGenNode> UnderGenPipeline::add_node_of_type(
    const String &class_name, const String &node_name,
    const Dictionary &properties, const Vector2 &editor_pos)
{
    Ref<UnderGenNode> node = ClassDB::instantiate(class_name);
    if (node.is_null()) {
        UtilityFunctions::printerr("UnderGenPipeline: Could not instantiate class '", class_name, "'.");
        return node;
    }

    node->set_name(node_name);
    node->set_editor_position(editor_pos);

    // Apply properties
    Array keys = properties.keys();
    for (int i = 0; i < keys.size(); i++) {
        String key = keys[i];
        Variant value = properties[key];
        node->set(key, value);
    }

    add_pipeline_node(node);
    return node;
}

void UnderGenPipeline::connect(const String &from_name, int from_port,
                                const String &to_name, int to_port)
{
    add_connection(from_name, from_port, to_name, to_port);
}

// ── Dictionary Serialization ─────────────────────────────────────────────────

Dictionary UnderGenPipeline::to_dictionary() const {
    Dictionary d;

    // Serialize nodes
    Array node_dicts;
    for (int i = 0; i < nodes.size(); i++) {
        Ref<UnderGenNode> n = nodes[i];
        if (n.is_null()) continue;

        Dictionary nd;
        nd["name"] = n->get_name();
        nd["type"] = n->get_class();
        nd["pos"]  = n->get_editor_position();

        // Collect all exported properties (that differ from defaults would be ideal,
        // but for clarity we collect all script-visible properties)
        Array prop_list = n->get_property_list();
        Dictionary node_props;
        for (int j = 0; j < prop_list.size(); j++) {
            Dictionary pinfo = prop_list[j];
            String pname = pinfo["name"];
            // Skip internal / metadata properties
            if (pname.begins_with("_") || pname == "name" || pname == "editor_position"
                || pname == "resource_local_to_scene" || pname == "script"
                || pname == "resource_path" || pname == "resource_name"
                || pname == "vox_spawn_entries" || pname == "vox_material_entries")
                continue;

            uint32_t usage = pinfo.get("usage", 0);
            // Only serialize properties marked as STORAGE
            if (!(usage & PROPERTY_USAGE_STORAGE)) continue;

            Variant val = n->get(pname);
            // Skip default values to keep the spec clean
            // (We'll include all for simplicity; users can trim.)
            node_props[pname] = val;
        }
        nd["properties"] = node_props;
        node_dicts.append(nd);
    }
    d["nodes"] = node_dicts;

    // Serialize connections
    d["connections"] = connections.duplicate();

    return d;
}

void UnderGenPipeline::from_dictionary(const Dictionary &d) {
    // Clear existing
    nodes.clear();
    connections.clear();

    // Rebuild nodes
    if (d.has("nodes")) {
        Array node_dicts = d["nodes"];
        for (int i = 0; i < node_dicts.size(); i++) {
            Dictionary nd = node_dicts[i];
            String type = nd.get("type", "");
            String name = nd.get("name", "");
            Vector2 pos  = nd.get("pos", Vector2());
            Dictionary props = nd.get("properties", Dictionary());

            Ref<UnderGenNode> node = add_node_of_type(type, name, props, pos);
            if (node.is_null()) {
                UtilityFunctions::printerr("UnderGenPipeline: Skipped unknown node type '", type, "'.");
            }
        }
    }

    // Rebuild connections
    if (d.has("connections")) {
        connections = d["connections"];
    }
}

// ── JSON Serialization ───────────────────────────────────────────────────────

String UnderGenPipeline::to_json(bool indent) const {
    Dictionary d = to_dictionary();
    return JSON::stringify(d, indent ? "    " : "");
}

Error UnderGenPipeline::from_json(const String &json_str,
                                    int *r_err_line, String *r_err_message) {
    Ref<JSON> json;
    json.instantiate();
    Error err = json->parse(json_str);
    if (err != OK) {
        if (r_err_line)    *r_err_line = json->get_error_line();
        if (r_err_message) *r_err_message = json->get_error_message();
        return err;
    }

    Variant parsed = json->get_data();
    if (parsed.get_type() != Variant::DICTIONARY) {
        if (r_err_line)    *r_err_line = 0;
        if (r_err_message) *r_err_message = "JSON root must be an object.";
        return ERR_PARSE_ERROR;
    }

    from_dictionary(parsed);
    return OK;
}

Error UnderGenPipeline::_from_json_bound(const String &json_str) {
    return from_json(json_str);
}

Error UnderGenPipeline::save_to_file(const String &path) {
    return ResourceSaver::get_singleton()->save(Ref<Resource>(this), path);
}

Ref<UnderGenPipeline> UnderGenPipeline::load_from_json_file(const String &path) {
    Ref<FileAccess> f = FileAccess::open(path, FileAccess::READ);
    if (f.is_null()) {
        UtilityFunctions::printerr("UnderGenPipeline: Could not open file: ", path);
        return Ref<UnderGenPipeline>();
    }
    String content = f->get_as_text();
    f->close();

    Ref<UnderGenPipeline> pipeline;
    pipeline.instantiate();
    Error err = pipeline->from_json(content);
    if (err != OK) {
        UtilityFunctions::printerr("UnderGenPipeline: JSON parse error in ", path, ": ", err);
        return Ref<UnderGenPipeline>();
    }
    return pipeline;
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
