// src/undergen_grammar_node.cpp
//
// Grammar expansion engine — a C++ port of grammar_generator.gd.
// Reads a LevelGrammarResource saved on disk via Object::get() so it
// works with GDScript @export properties without needing a C++ resource class.
//
// Algorithm (matches grammar_generator.gd):
//   1. Start with a single "root" node whose symbol == grammar.axiom.
//   2. For each iteration: for every node in the graph, look up matching rules,
//      pick one by weighted probability, and replace the node with the rule's
//      RHS subgraph.  Rewire external edges through the replacement's entry/exit.
//   3. Output: Dictionary{ "nodes": Array, "edges": Array } on port 0.

#include "undergen_grammar_node.h"

#include <godot_cpp/classes/resource_loader.hpp>
#include <godot_cpp/variant/utility_functions.hpp>

#include <vector>
#include <map>
#include <cmath>

namespace godot {

// ─── constructor ────────────────────────────────────────────────────────────

UnderGenGrammarNode::UnderGenGrammarNode() {
    rng.instantiate();
}

// ─── bindings ───────────────────────────────────────────────────────────────

void UnderGenGrammarNode::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_grammar_resource_path", "path"),
                         &UnderGenGrammarNode::set_grammar_resource_path);
    ClassDB::bind_method(D_METHOD("get_grammar_resource_path"),
                         &UnderGenGrammarNode::get_grammar_resource_path);

    ClassDB::bind_method(D_METHOD("set_iterations", "iterations"),
                         &UnderGenGrammarNode::set_iterations);
    ClassDB::bind_method(D_METHOD("get_iterations"),
                         &UnderGenGrammarNode::get_iterations);

    ClassDB::bind_method(D_METHOD("set_max_nodes", "max"),
                         &UnderGenGrammarNode::set_max_nodes);
    ClassDB::bind_method(D_METHOD("get_max_nodes"),
                         &UnderGenGrammarNode::get_max_nodes);

    ADD_PROPERTY(
        PropertyInfo(Variant::STRING, "grammar_resource_path",
                     PROPERTY_HINT_FILE, "*.tres,*.res"),
        "set_grammar_resource_path", "get_grammar_resource_path");

    ADD_PROPERTY(
        PropertyInfo(Variant::INT, "iterations",
                     PROPERTY_HINT_RANGE, "1,32,1"),
        "set_iterations", "get_iterations");

    ADD_PROPERTY(
        PropertyInfo(Variant::INT, "max_nodes",
                     PROPERTY_HINT_RANGE, "1,1000,1"),
        "set_max_nodes", "get_max_nodes");
}

// ─── accessors ──────────────────────────────────────────────────────────────

void   UnderGenGrammarNode::set_grammar_resource_path(const String& p) { grammar_resource_path = p; }
String UnderGenGrammarNode::get_grammar_resource_path() const          { return grammar_resource_path; }
void   UnderGenGrammarNode::set_iterations(int v) { iterations = v; }
int    UnderGenGrammarNode::get_iterations() const { return iterations; }
void   UnderGenGrammarNode::set_max_nodes(int v)  { max_nodes = v; }
int    UnderGenGrammarNode::get_max_nodes() const  { return max_nodes; }

// ─── internal helpers ───────────────────────────────────────────────────────

bool UnderGenGrammarNode::_eval_condition(
    const String& var, const String& op, double val,
    const Dictionary& state) const
{
    if (var.is_empty()) return true; // no condition = always fire
    double cur = (double)(state.get(var, 0.0));
    if (op == "<")  return cur < val;
    if (op == ">")  return cur > val;
    if (op == "<=") return cur <= val;
    if (op == ">=") return cur >= val;
    if (op == "==") return std::fabs(cur - val) < 1e-6;
    if (op == "!=") return std::fabs(cur - val) >= 1e-6;
    return true;
}

Ref<Resource> UnderGenGrammarNode::_pick_rule(
    const Array& candidates, const Dictionary& state)
{
    std::vector<Ref<Resource>> valid;
    for (int i = 0; i < candidates.size(); ++i) {
        Ref<Resource> rule = candidates[i];
        if (!rule.is_valid()) continue;
        String cvar = rule->get("condition_var");
        String cop  = rule->get("condition_op");
        double cval = (double)(rule->get("condition_val"));
        if (_eval_condition(cvar, cop, cval, state))
            valid.push_back(rule);
    }
    if (valid.empty()) return Ref<Resource>();

    // Weighted random pick
    double total = 0.0;
    for (auto& r : valid) total += (double)(r->get("probability"));
    if (total <= 0.0) return valid[0];

    double pick = (double)rng->randf() * total;
    double accum = 0.0;
    for (auto& r : valid) {
        accum += (double)(r->get("probability"));
        if (pick <= accum) return r;
    }
    return valid.back();
}

Dictionary UnderGenGrammarNode::_make_node(
    const String& id, const String& symbol,
    const Dictionary& constraints,
    const std::map<String, Dictionary>& rt_map) const
{
    Dictionary n;
    n["id"]          = id;
    n["symbol"]      = symbol;
    n["type"]        = symbol.to_lower(); // zone label for carver
    n["constraints"] = constraints;

    auto it = rt_map.find(symbol);
    if (it != rt_map.end()) {
        const Dictionary& rt = it->second;
        n["min_size"] = rt.get("min_size", Vector3i(5, 3, 5));
        n["max_size"] = rt.get("max_size", Vector3i(10, 6, 10));
        String vox = rt.get("vox_path", "");
        if (!vox.is_empty()) n["vox_path"] = vox;
        n["exclude_from_smoothing"] = rt.get("exclude_from_smoothing", false);
        n["exclude_from_warping"] = rt.get("exclude_from_warping", false);
    } else {
        n["min_size"] = Vector3i(5, 3, 5);
        n["max_size"] = Vector3i(10, 6, 10);
    }
    return n;
}

// ─── main execute ───────────────────────────────────────────────────────────

void UnderGenGrammarNode::_execute(const Dictionary& inputs, Dictionary& outputs) {
    if (grammar_resource_path.is_empty()) {
        UtilityFunctions::printerr("UnderGenGrammarNode: grammar_resource_path is empty.");
        return;
    }

    // Load grammar resource (GDScript LevelGrammarResource saved as .tres)
    Ref<Resource> grammar = ResourceLoader::get_singleton()->load(grammar_resource_path);
    if (!grammar.is_valid()) {
        UtilityFunctions::printerr(
            "UnderGenGrammarNode: Could not load grammar: ", grammar_resource_path);
        return;
    }

    // Seed the RNG from port 0
    int64_t seed = inputs.get(0, (int64_t)12345);
    rng->set_seed(seed);

    // Read grammar properties via Godot object API
    String axiom = grammar->get("axiom");

    // Build room-type lookup map: symbol -> { min_size, max_size, vox_path, ... }
    std::map<String, Dictionary> rt_map;
    Array room_types_arr = grammar->get("room_types");
    for (int i = 0; i < room_types_arr.size(); ++i) {
        Ref<Resource> rt = room_types_arr[i];
        if (!rt.is_valid()) continue;
        String sym = rt->get("symbol");
        Dictionary d;
        d["min_size"] = rt->get("min_size");
        d["max_size"] = rt->get("max_size");
        d["vox_path"] = rt->get("vox_path");
        d["exclude_from_smoothing"] = rt->get("exclude_from_smoothing");
        d["exclude_from_warping"] = rt->get("exclude_from_warping");
        rt_map[sym] = d;
    }

    // Build rules index: lhs_symbol -> [rule, rule, ...]
    std::map<String, Array> rules_by_sym;
    Array rules_arr = grammar->get("rules");
    for (int i = 0; i < rules_arr.size(); ++i) {
        Ref<Resource> rule = rules_arr[i];
        if (!rule.is_valid()) continue;
        String lhs = rule->get("lhs_symbol");
        rules_by_sym[lhs].append(rule);
    }

    // Initial state (can be seeded from outside later)
    Dictionary state;

    // ── Seed graph: one root node with the axiom symbol ──────────────────────
    Array nodes, edges;
    Dictionary root_c;
    root_c["fixed_pos"] = Vector3(0.0f, 0.0f, 0.0f);
    nodes.append(_make_node("root", axiom, root_c, rt_map));

    // ── Iterative expansion ──────────────────────────────────────────────────
    for (int iter = 0; iter < iterations; ++iter) {
        if ((int)nodes.size() >= max_nodes) {
            UtilityFunctions::print(
                "UnderGenGrammarNode: max_nodes (", max_nodes, ") reached, stopping.");
            break;
        }

        std::map<String, Dictionary> replacements; // old_id -> {entry,exit,nodes,edges}
        Array preserved;
        bool any_change = false;

        for (int ni = 0; ni < nodes.size(); ++ni) {
            Dictionary node = nodes[ni];
            String symbol   = node.get("symbol", "");
            String node_id  = node.get("id", "");

            auto it = rules_by_sym.find(symbol);
            if (it == rules_by_sym.end()) {
                preserved.append(node);
                continue;
            }

            Ref<Resource> rule = _pick_rule(it->second, state);
            if (!rule.is_valid()) {
                preserved.append(node);
                continue;
            }

            any_change = true;

            // ── Apply state actions ──────────────────────────────────────────
            Array actions = rule->get("actions");
            for (int ai = 0; ai < actions.size(); ++ai) {
                Dictionary act = actions[ai];
                String var = act.get("var", "");
                double delta = (double)(act.get("delta", 0.0));
                if (!var.is_empty()) {
                    double cur = (double)(state.get(var, 0.0));
                    state[var] = cur + delta;
                }
            }

            // ── Build replacement subgraph ───────────────────────────────────
            Array sub_nodes, sub_edges;
            std::map<String, String> id_map; // local_id -> global_unique_id

            Array rhs_nodes_arr = rule->get("rhs_nodes");
            for (int rni = 0; rni < rhs_nodes_arr.size(); ++rni) {
                Dictionary rhs_n   = rhs_nodes_arr[rni];
                String local_id    = rhs_n.get("id", "");
                String unique_id   = node_id + "-" + local_id + "-" + String::num_int64(iter);
                id_map[local_id]   = unique_id;

                String sym         = rhs_n.get("symbol", "");
                Dictionary cons    = ((Dictionary)rhs_n.get("constraints", Dictionary())).duplicate();

                // Resolve internal relative_to references
                if (cons.has("relative_to")) {
                    String rel = cons.get("relative_to", "");
                    auto rid = id_map.find(rel);
                    if (rid != id_map.end()) cons["relative_to"] = rid->second;
                }

                sub_nodes.append(_make_node(unique_id, sym, cons, rt_map));
            }

            Array rhs_edges_arr = rule->get("rhs_edges");
            for (int rei = 0; rei < rhs_edges_arr.size(); ++rei) {
                Dictionary rhs_e = rhs_edges_arr[rei];
                String fl = rhs_e.get("from", "");
                String tl = rhs_e.get("to",   "");
                auto fi = id_map.find(fl), ti = id_map.find(tl);
                if (fi != id_map.end() && ti != id_map.end()) {
                    Dictionary e;
                    e["from"] = fi->second;
                    e["to"]   = ti->second;
                    e["type"] = rhs_e.get("type", "corridor");
                    sub_edges.append(e);
                }
            }

            // ── Determine entry/exit in new subgraph ─────────────────────────
            String entry_global, exit_global;

            String entry_local = rule->get("entry_node_id");
            String exit_local  = rule->get("exit_node_id");

            if (!entry_local.is_empty() && id_map.count(entry_local)) {
                entry_global = id_map[entry_local];
            } else if (sub_nodes.size() > 0) {
                entry_global = ((Dictionary)sub_nodes[0]).get("id", node_id);
            } else {
                entry_global = node_id;
            }

            if (!exit_local.is_empty() && id_map.count(exit_local)) {
                exit_global = id_map[exit_local];
            } else {
                exit_global = entry_global;
            }

            Dictionary rep;
            rep["entry"] = entry_global;
            rep["exit"]  = exit_global;
            rep["nodes"] = sub_nodes;
            rep["edges"] = sub_edges;
            replacements[node_id] = rep;
        }

        if (!any_change) break;

        // ── Rebuild graph ────────────────────────────────────────────────────
        Array next_nodes, next_edges;
        next_nodes.append_array(preserved);

        for (auto& [id, rep] : replacements) {
            next_nodes.append_array(rep.get("nodes", Array()));
            next_edges.append_array(rep.get("edges", Array()));
        }

        // Rewire edges through entry/exit of replaced nodes
        for (int ei = 0; ei < edges.size(); ++ei) {
            Dictionary old_e = edges[ei];
            String from = old_e.get("from", "");
            String to   = old_e.get("to",   "");

            String new_from = from;
            auto fit = replacements.find(from);
            if (fit != replacements.end())
                new_from = fit->second.get("exit", from);

            String new_to = to;
            auto tit = replacements.find(to);
            if (tit != replacements.end())
                new_to = tit->second.get("entry", to);

            if (new_from != new_to) {
                Dictionary ne;
                ne["from"] = new_from;
                ne["to"]   = new_to;
                ne["type"] = old_e.get("type", "corridor");
                next_edges.append(ne);
            }
        }

        nodes = next_nodes;
        edges = next_edges;
    }

    // ── Output ───────────────────────────────────────────────────────────────
    Dictionary result;
    result["nodes"] = nodes;
    result["edges"] = edges;
    outputs[0] = result;

    UtilityFunctions::print(
        "UnderGenGrammarNode: expanded to ", nodes.size(),
        " nodes, ", edges.size(), " edges.");
}

Dictionary UnderGenGrammarNode::get_pipeline_input_defaults(const Dictionary &global_inputs) const {
    Dictionary defaults;
    if (global_inputs.has(0)) defaults[0] = global_inputs[0];
    return defaults;
}

} // namespace godot
