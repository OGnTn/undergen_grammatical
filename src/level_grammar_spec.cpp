// src/level_grammar_spec.cpp
//
// C++ implementation of the grammar specification resources.
// These mirror the GDScript LevelGrammarResource / GraphRule / RoomType
// property-by-property so they are drop-in compatible with the existing
// UnderGenGrammarNode pipeline node.

#include "level_grammar_spec.h"

#include <godot_cpp/classes/json.hpp>
#include <godot_cpp/classes/file_access.hpp>
#include <godot_cpp/classes/resource_saver.hpp>
#include <godot_cpp/classes/resource_loader.hpp>
#include <godot_cpp/variant/utility_functions.hpp>

namespace godot {

// ═══════════════════════════════════════════════════════════════════════════
//  LevelGrammarRoomTypeSpec
// ═══════════════════════════════════════════════════════════════════════════

LevelGrammarRoomTypeSpec::LevelGrammarRoomTypeSpec() = default;

void LevelGrammarRoomTypeSpec::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_symbol", "symbol"), &LevelGrammarRoomTypeSpec::set_symbol);
    ClassDB::bind_method(D_METHOD("get_symbol"), &LevelGrammarRoomTypeSpec::get_symbol);
    ClassDB::bind_method(D_METHOD("set_color", "color"), &LevelGrammarRoomTypeSpec::set_color);
    ClassDB::bind_method(D_METHOD("get_color"), &LevelGrammarRoomTypeSpec::get_color);
    ClassDB::bind_method(D_METHOD("set_weight", "weight"), &LevelGrammarRoomTypeSpec::set_weight);
    ClassDB::bind_method(D_METHOD("get_weight"), &LevelGrammarRoomTypeSpec::get_weight);
    ClassDB::bind_method(D_METHOD("set_min_size", "min_size"), &LevelGrammarRoomTypeSpec::set_min_size);
    ClassDB::bind_method(D_METHOD("get_min_size"), &LevelGrammarRoomTypeSpec::get_min_size);
    ClassDB::bind_method(D_METHOD("set_max_size", "max_size"), &LevelGrammarRoomTypeSpec::set_max_size);
    ClassDB::bind_method(D_METHOD("get_max_size"), &LevelGrammarRoomTypeSpec::get_max_size);
    ClassDB::bind_method(D_METHOD("set_vox_path", "vox_path"), &LevelGrammarRoomTypeSpec::set_vox_path);
    ClassDB::bind_method(D_METHOD("get_vox_path"), &LevelGrammarRoomTypeSpec::get_vox_path);
    ClassDB::bind_method(D_METHOD("set_exclude_from_smoothing", "exclude"), &LevelGrammarRoomTypeSpec::set_exclude_from_smoothing);
    ClassDB::bind_method(D_METHOD("get_exclude_from_smoothing"), &LevelGrammarRoomTypeSpec::get_exclude_from_smoothing);
    ClassDB::bind_method(D_METHOD("set_exclude_from_warping", "exclude"), &LevelGrammarRoomTypeSpec::set_exclude_from_warping);
    ClassDB::bind_method(D_METHOD("get_exclude_from_warping"), &LevelGrammarRoomTypeSpec::get_exclude_from_warping);

    ClassDB::bind_method(D_METHOD("to_dictionary"), &LevelGrammarRoomTypeSpec::to_dictionary);
    ClassDB::bind_method(D_METHOD("from_dictionary", "d"), &LevelGrammarRoomTypeSpec::from_dictionary);

    ADD_PROPERTY(PropertyInfo(Variant::STRING, "symbol"), "set_symbol", "get_symbol");
    ADD_PROPERTY(PropertyInfo(Variant::COLOR, "color"), "set_color", "get_color");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "weight"), "set_weight", "get_weight");
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3I, "min_size"), "set_min_size", "get_min_size");
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3I, "max_size"), "set_max_size", "get_max_size");
    ADD_PROPERTY(PropertyInfo(Variant::STRING, "vox_path", PROPERTY_HINT_FILE, "*.vox,*.tres,*.res"), "set_vox_path", "get_vox_path");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "exclude_from_smoothing"), "set_exclude_from_smoothing", "get_exclude_from_smoothing");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "exclude_from_warping"), "set_exclude_from_warping", "get_exclude_from_warping");
}

// ── Accessors ─────────────────────────────────────────────────────────────

void LevelGrammarRoomTypeSpec::set_symbol(const String &p_sym) { symbol = p_sym; }
String LevelGrammarRoomTypeSpec::get_symbol() const { return symbol; }
void LevelGrammarRoomTypeSpec::set_color(const Color &p_color) { color = p_color; }
Color LevelGrammarRoomTypeSpec::get_color() const { return color; }
void LevelGrammarRoomTypeSpec::set_weight(float p_w) { weight = p_w; }
float LevelGrammarRoomTypeSpec::get_weight() const { return weight; }
void LevelGrammarRoomTypeSpec::set_min_size(const Vector3i &p_size) { min_size = p_size; }
Vector3i LevelGrammarRoomTypeSpec::get_min_size() const { return min_size; }
void LevelGrammarRoomTypeSpec::set_max_size(const Vector3i &p_size) { max_size = p_size; }
Vector3i LevelGrammarRoomTypeSpec::get_max_size() const { return max_size; }
void LevelGrammarRoomTypeSpec::set_vox_path(const String &p_path) { vox_path = p_path; }
String LevelGrammarRoomTypeSpec::get_vox_path() const { return vox_path; }
void LevelGrammarRoomTypeSpec::set_exclude_from_smoothing(bool p_exclude) { exclude_from_smoothing = p_exclude; }
bool LevelGrammarRoomTypeSpec::get_exclude_from_smoothing() const { return exclude_from_smoothing; }
void LevelGrammarRoomTypeSpec::set_exclude_from_warping(bool p_exclude) { exclude_from_warping = p_exclude; }
bool LevelGrammarRoomTypeSpec::get_exclude_from_warping() const { return exclude_from_warping; }

// ── Helper: safely extract Vector3i from a Variant ───────────────────────
// JSON arrays like [4,3,4] come through as Godot Array, not Vector3i.
static Vector3i _variant_to_vector3i(const Variant &v, const Vector3i &fallback = Vector3i()) {
    if (v.get_type() == Variant::VECTOR3I) {
        return Vector3i(v);
    }
    if (v.get_type() == Variant::VECTOR3) {
        Vector3 fv = v;
        return Vector3i((int)fv.x, (int)fv.y, (int)fv.z);
    }
    if (v.get_type() == Variant::ARRAY) {
        Array arr = v;
        if (arr.size() >= 3) {
            return Vector3i(
                (int)(double)(arr[0]),
                (int)(double)(arr[1]),
                (int)(double)(arr[2])
            );
        }
    }
    UtilityFunctions::printerr("LevelGrammarSpec: Cannot convert Variant (type=", v.get_type(), ") to Vector3i, using fallback ", fallback);
    return fallback;
}

// ── Serialization ─────────────────────────────────────────────────────────

Dictionary LevelGrammarRoomTypeSpec::to_dictionary() const {
    Dictionary d;
    d["symbol"]   = symbol;
    d["color"]    = color;
    d["weight"]   = weight;
    d["min_size"] = min_size;
    d["max_size"] = max_size;
    if (!vox_path.is_empty())
        d["vox_path"] = vox_path;
    if (exclude_from_smoothing)
        d["exclude_from_smoothing"] = exclude_from_smoothing;
    if (exclude_from_warping)
        d["exclude_from_warping"] = exclude_from_warping;
    return d;
}

void LevelGrammarRoomTypeSpec::from_dictionary(const Dictionary &d) {
    if (d.has("symbol"))   set_symbol(d["symbol"]);
    if (d.has("color"))    set_color(d["color"]);
    if (d.has("weight"))   set_weight(d["weight"]);
    if (d.has("min_size")) set_min_size(_variant_to_vector3i(d["min_size"], Vector3i(5, 3, 5)));
    if (d.has("max_size")) set_max_size(_variant_to_vector3i(d["max_size"], Vector3i(10, 6, 10)));
    if (d.has("vox_path")) set_vox_path(d["vox_path"]);
    if (d.has("exclude_from_smoothing")) set_exclude_from_smoothing(d["exclude_from_smoothing"]);
    if (d.has("exclude_from_warping")) set_exclude_from_warping(d["exclude_from_warping"]);
}


// ═══════════════════════════════════════════════════════════════════════════
//  LevelGrammarRuleSpec
// ═══════════════════════════════════════════════════════════════════════════

LevelGrammarRuleSpec::LevelGrammarRuleSpec() = default;

void LevelGrammarRuleSpec::_bind_methods() {
    // Core
    ClassDB::bind_method(D_METHOD("set_rule_name", "name"), &LevelGrammarRuleSpec::set_rule_name);
    ClassDB::bind_method(D_METHOD("get_rule_name"), &LevelGrammarRuleSpec::get_rule_name);
    ClassDB::bind_method(D_METHOD("set_lhs_symbol", "symbol"), &LevelGrammarRuleSpec::set_lhs_symbol);
    ClassDB::bind_method(D_METHOD("get_lhs_symbol"), &LevelGrammarRuleSpec::get_lhs_symbol);
    ClassDB::bind_method(D_METHOD("set_probability", "prob"), &LevelGrammarRuleSpec::set_probability);
    ClassDB::bind_method(D_METHOD("get_probability"), &LevelGrammarRuleSpec::get_probability);
    ClassDB::bind_method(D_METHOD("set_entry_node_id", "id"), &LevelGrammarRuleSpec::set_entry_node_id);
    ClassDB::bind_method(D_METHOD("get_entry_node_id"), &LevelGrammarRuleSpec::get_entry_node_id);
    ClassDB::bind_method(D_METHOD("set_exit_node_id", "id"), &LevelGrammarRuleSpec::set_exit_node_id);
    ClassDB::bind_method(D_METHOD("get_exit_node_id"), &LevelGrammarRuleSpec::get_exit_node_id);

    // Condition
    ClassDB::bind_method(D_METHOD("set_condition_var", "var"), &LevelGrammarRuleSpec::set_condition_var);
    ClassDB::bind_method(D_METHOD("get_condition_var"), &LevelGrammarRuleSpec::get_condition_var);
    ClassDB::bind_method(D_METHOD("set_condition_op", "op"), &LevelGrammarRuleSpec::set_condition_op);
    ClassDB::bind_method(D_METHOD("get_condition_op"), &LevelGrammarRuleSpec::get_condition_op);
    ClassDB::bind_method(D_METHOD("set_condition_val", "val"), &LevelGrammarRuleSpec::set_condition_val);
    ClassDB::bind_method(D_METHOD("get_condition_val"), &LevelGrammarRuleSpec::get_condition_val);

    // Actions & RHS
    ClassDB::bind_method(D_METHOD("set_actions", "actions"), &LevelGrammarRuleSpec::set_actions);
    ClassDB::bind_method(D_METHOD("get_actions"), &LevelGrammarRuleSpec::get_actions);
    ClassDB::bind_method(D_METHOD("set_rhs_nodes", "nodes"), &LevelGrammarRuleSpec::set_rhs_nodes);
    ClassDB::bind_method(D_METHOD("get_rhs_nodes"), &LevelGrammarRuleSpec::get_rhs_nodes);
    ClassDB::bind_method(D_METHOD("set_rhs_edges", "edges"), &LevelGrammarRuleSpec::set_rhs_edges);
    ClassDB::bind_method(D_METHOD("get_rhs_edges"), &LevelGrammarRuleSpec::get_rhs_edges);

    ClassDB::bind_method(D_METHOD("to_dictionary"), &LevelGrammarRuleSpec::to_dictionary);
    ClassDB::bind_method(D_METHOD("from_dictionary", "d"), &LevelGrammarRuleSpec::from_dictionary);

    ADD_PROPERTY(PropertyInfo(Variant::STRING, "rule_name"), "set_rule_name", "get_rule_name");
    ADD_PROPERTY(PropertyInfo(Variant::STRING, "lhs_symbol"), "set_lhs_symbol", "get_lhs_symbol");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "probability", PROPERTY_HINT_RANGE, "0.0,100.0,0.1"), "set_probability", "get_probability");
    ADD_PROPERTY(PropertyInfo(Variant::STRING, "entry_node_id"), "set_entry_node_id", "get_entry_node_id");
    ADD_PROPERTY(PropertyInfo(Variant::STRING, "exit_node_id"), "set_exit_node_id", "get_exit_node_id");

    ADD_GROUP("Condition", "condition_");
    ADD_PROPERTY(PropertyInfo(Variant::STRING, "condition_var"), "set_condition_var", "get_condition_var");
    ADD_PROPERTY(PropertyInfo(Variant::STRING, "condition_op"), "set_condition_op", "get_condition_op");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "condition_val"), "set_condition_val", "get_condition_val");

    ADD_GROUP("Actions", "");
    ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "actions"), "set_actions", "get_actions");

    ADD_GROUP("Graph Data", "");
    ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "rhs_nodes"), "set_rhs_nodes", "get_rhs_nodes");
    ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "rhs_edges"), "set_rhs_edges", "get_rhs_edges");
}

// ── Accessors ─────────────────────────────────────────────────────────────

void LevelGrammarRuleSpec::set_rule_name(const String &p) { rule_name = p; }
String LevelGrammarRuleSpec::get_rule_name() const { return rule_name; }
void LevelGrammarRuleSpec::set_lhs_symbol(const String &p) { lhs_symbol = p; }
String LevelGrammarRuleSpec::get_lhs_symbol() const { return lhs_symbol; }
void LevelGrammarRuleSpec::set_probability(float p) { probability = p; }
float LevelGrammarRuleSpec::get_probability() const { return probability; }
void LevelGrammarRuleSpec::set_entry_node_id(const String &p) { entry_node_id = p; }
String LevelGrammarRuleSpec::get_entry_node_id() const { return entry_node_id; }
void LevelGrammarRuleSpec::set_exit_node_id(const String &p) { exit_node_id = p; }
String LevelGrammarRuleSpec::get_exit_node_id() const { return exit_node_id; }
void LevelGrammarRuleSpec::set_condition_var(const String &p) { condition_var = p; }
String LevelGrammarRuleSpec::get_condition_var() const { return condition_var; }
void LevelGrammarRuleSpec::set_condition_op(const String &p) { condition_op = p; }
String LevelGrammarRuleSpec::get_condition_op() const { return condition_op; }
void LevelGrammarRuleSpec::set_condition_val(float p) { condition_val = p; }
float LevelGrammarRuleSpec::get_condition_val() const { return condition_val; }
void LevelGrammarRuleSpec::set_actions(const Array &p) { actions = p; }
Array LevelGrammarRuleSpec::get_actions() const { return actions; }
void LevelGrammarRuleSpec::set_rhs_nodes(const Array &p) { rhs_nodes = p; }
Array LevelGrammarRuleSpec::get_rhs_nodes() const { return rhs_nodes; }
void LevelGrammarRuleSpec::set_rhs_edges(const Array &p) { rhs_edges = p; }
Array LevelGrammarRuleSpec::get_rhs_edges() const { return rhs_edges; }

// ── Serialization ─────────────────────────────────────────────────────────

Dictionary LevelGrammarRuleSpec::to_dictionary() const {
    Dictionary d;
    d["rule_name"]    = rule_name;
    d["lhs_symbol"]   = lhs_symbol;
    d["probability"]  = probability;

    if (!entry_node_id.is_empty()) d["entry_node_id"] = entry_node_id;
    if (!exit_node_id.is_empty())  d["exit_node_id"]  = exit_node_id;

    // Condition (only emit if a var is set)
    if (!condition_var.is_empty()) {
        d["condition_var"] = condition_var;
        d["condition_op"]  = condition_op;
        d["condition_val"] = condition_val;
    }

    if (!actions.is_empty())
        d["actions"] = actions;

    d["rhs_nodes"] = rhs_nodes;
    d["rhs_edges"] = rhs_edges;

    return d;
}

void LevelGrammarRuleSpec::from_dictionary(const Dictionary &d) {
    if (d.has("rule_name"))    set_rule_name(d["rule_name"]);
    if (d.has("lhs_symbol"))   set_lhs_symbol(d["lhs_symbol"]);
    if (d.has("probability"))  set_probability(d["probability"]);
    if (d.has("entry_node_id")) set_entry_node_id(d["entry_node_id"]);
    if (d.has("exit_node_id"))  set_exit_node_id(d["exit_node_id"]);
    if (d.has("condition_var")) set_condition_var(d["condition_var"]);
    if (d.has("condition_op"))  set_condition_op(d["condition_op"]);
    if (d.has("condition_val")) set_condition_val(d["condition_val"]);
    if (d.has("actions"))      set_actions(d["actions"]);
    if (d.has("rhs_nodes"))    set_rhs_nodes(d["rhs_nodes"]);
    if (d.has("rhs_edges"))    set_rhs_edges(d["rhs_edges"]);
}


// ═══════════════════════════════════════════════════════════════════════════
//  LevelGrammarSpec
// ═══════════════════════════════════════════════════════════════════════════

LevelGrammarSpec::LevelGrammarSpec() = default;

void LevelGrammarSpec::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_axiom", "axiom"), &LevelGrammarSpec::set_axiom);
    ClassDB::bind_method(D_METHOD("get_axiom"), &LevelGrammarSpec::get_axiom);
    ClassDB::bind_method(D_METHOD("set_state_variables", "vars"), &LevelGrammarSpec::set_state_variables);
    ClassDB::bind_method(D_METHOD("get_state_variables"), &LevelGrammarSpec::get_state_variables);
    ClassDB::bind_method(D_METHOD("set_room_types", "types"), &LevelGrammarSpec::set_room_types);
    ClassDB::bind_method(D_METHOD("get_room_types"), &LevelGrammarSpec::get_room_types);
    ClassDB::bind_method(D_METHOD("set_rules", "rules"), &LevelGrammarSpec::set_rules);
    ClassDB::bind_method(D_METHOD("get_rules"), &LevelGrammarSpec::get_rules);

    // Builder
    ClassDB::bind_method(D_METHOD("add_room_type", "symbol", "color", "weight", "min_size", "max_size", "vox_path"),
                         &LevelGrammarSpec::add_room_type,
                         DEFVAL(Color(0.4f, 0.6f, 0.9f)), DEFVAL(1.0f),
                         DEFVAL(Vector3i(5, 3, 5)), DEFVAL(Vector3i(10, 6, 10)), DEFVAL(""));
    ClassDB::bind_method(D_METHOD("add_room_type_vox", "symbol", "vox_path", "color"),
                         &LevelGrammarSpec::add_room_type_vox,
                         DEFVAL(Color(0.4f, 0.6f, 0.9f)));
    ClassDB::bind_method(D_METHOD("add_rule", "rule_name", "lhs_symbol", "probability", "rhs_nodes", "rhs_edges", "entry_node_id", "exit_node_id"),
                         &LevelGrammarSpec::add_rule,
                         DEFVAL(""), DEFVAL(""));
    ClassDB::bind_method(D_METHOD("add_rule_conditional", "rule_name", "lhs_symbol", "probability", "condition_var", "condition_op", "condition_val", "actions", "rhs_nodes", "rhs_edges", "entry_node_id", "exit_node_id"),
                         &LevelGrammarSpec::add_rule_conditional,
                         DEFVAL(""), DEFVAL(""));

    // Serialization
    ClassDB::bind_method(D_METHOD("to_dictionary"), &LevelGrammarSpec::to_dictionary);
    ClassDB::bind_method(D_METHOD("from_dictionary", "d"), &LevelGrammarSpec::from_dictionary);
    ClassDB::bind_method(D_METHOD("to_json", "indent"), &LevelGrammarSpec::to_json, DEFVAL(true));
    ClassDB::bind_method(D_METHOD("from_json", "json_str"), &LevelGrammarSpec::_from_json_bound);
    ClassDB::bind_method(D_METHOD("save_to_file", "path"), &LevelGrammarSpec::save_to_file);
    ClassDB::bind_static_method("LevelGrammarSpec", D_METHOD("load_from_json_file", "path"), &LevelGrammarSpec::load_from_json_file);

    ADD_PROPERTY(PropertyInfo(Variant::STRING, "axiom"), "set_axiom", "get_axiom");
    ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "state_variables"), "set_state_variables", "get_state_variables");
    ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "room_types"), "set_room_types", "get_room_types");
    ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "rules"), "set_rules", "get_rules");
}

// ── Accessors ─────────────────────────────────────────────────────────────

void LevelGrammarSpec::set_axiom(const String &p) { axiom = p; }
String LevelGrammarSpec::get_axiom() const { return axiom; }
void LevelGrammarSpec::set_state_variables(const Array &p) { state_variables = p; }
Array LevelGrammarSpec::get_state_variables() const { return state_variables; }
void LevelGrammarSpec::set_room_types(const TypedArray<LevelGrammarRoomTypeSpec> &p) { room_types = p; }
TypedArray<LevelGrammarRoomTypeSpec> LevelGrammarSpec::get_room_types() const { return room_types; }
void LevelGrammarSpec::set_rules(const TypedArray<LevelGrammarRuleSpec> &p) { rules = p; }
TypedArray<LevelGrammarRuleSpec> LevelGrammarSpec::get_rules() const { return rules; }

// ── Builder API ───────────────────────────────────────────────────────────

Ref<LevelGrammarRoomTypeSpec> LevelGrammarSpec::add_room_type(
    const String &symbol, const Color &color, float weight,
    const Vector3i &min_size, const Vector3i &max_size, const String &vox_path)
{
    Ref<LevelGrammarRoomTypeSpec> rt;
    rt.instantiate();
    rt->set_symbol(symbol);
    rt->set_color(color);
    rt->set_weight(weight);
    rt->set_min_size(min_size);
    rt->set_max_size(max_size);
    if (!vox_path.is_empty()) rt->set_vox_path(vox_path);
    room_types.append(rt);
    return rt;
}

Ref<LevelGrammarRoomTypeSpec> LevelGrammarSpec::add_room_type_vox(
    const String &symbol, const String &vox_path, const Color &color)
{
    return add_room_type(symbol, color, 1.0f,
                         Vector3i(5, 3, 5), Vector3i(10, 6, 10), vox_path);
}

Ref<LevelGrammarRuleSpec> LevelGrammarSpec::add_rule(
    const String &rule_name, const String &lhs_symbol, float probability,
    const Array &rhs_nodes, const Array &rhs_edges,
    const String &entry_node_id, const String &exit_node_id)
{
    Ref<LevelGrammarRuleSpec> r;
    r.instantiate();
    r->set_rule_name(rule_name);
    r->set_lhs_symbol(lhs_symbol);
    r->set_probability(probability);
    r->set_rhs_nodes(rhs_nodes);
    r->set_rhs_edges(rhs_edges);
    if (!entry_node_id.is_empty()) r->set_entry_node_id(entry_node_id);
    if (!exit_node_id.is_empty())  r->set_exit_node_id(exit_node_id);
    rules.append(r);
    return r;
}

Ref<LevelGrammarRuleSpec> LevelGrammarSpec::add_rule_conditional(
    const String &rule_name, const String &lhs_symbol, float probability,
    const String &condition_var, const String &condition_op, float condition_val,
    const Array &actions, const Array &rhs_nodes, const Array &rhs_edges,
    const String &entry_node_id, const String &exit_node_id)
{
    Ref<LevelGrammarRuleSpec> r;
    r.instantiate();
    r->set_rule_name(rule_name);
    r->set_lhs_symbol(lhs_symbol);
    r->set_probability(probability);
    r->set_condition_var(condition_var);
    r->set_condition_op(condition_op);
    r->set_condition_val(condition_val);
    r->set_actions(actions);
    r->set_rhs_nodes(rhs_nodes);
    r->set_rhs_edges(rhs_edges);
    if (!entry_node_id.is_empty()) r->set_entry_node_id(entry_node_id);
    if (!exit_node_id.is_empty())  r->set_exit_node_id(exit_node_id);
    rules.append(r);
    return r;
}

// ── Dictionary Serialization ──────────────────────────────────────────────

Dictionary LevelGrammarSpec::to_dictionary() const {
    Dictionary d;
    d["axiom"] = axiom;

    // State variables
    if (!state_variables.is_empty())
        d["state_variables"] = state_variables;

    // Room types
    Array rt_arr;
    for (int i = 0; i < room_types.size(); i++) {
        Ref<LevelGrammarRoomTypeSpec> rt = room_types[i];
        if (rt.is_valid())
            rt_arr.append(rt->to_dictionary());
    }
    d["room_types"] = rt_arr;

    // Rules
    Array r_arr;
    for (int i = 0; i < rules.size(); i++) {
        Ref<LevelGrammarRuleSpec> r = rules[i];
        if (r.is_valid())
            r_arr.append(r->to_dictionary());
    }
    d["rules"] = r_arr;

    return d;
}

void LevelGrammarSpec::from_dictionary(const Dictionary &d) {
    // Clear existing data
    room_types.clear();
    rules.clear();
    state_variables.clear();

    if (d.has("axiom")) set_axiom(d["axiom"]);

    if (d.has("state_variables")) {
        Array sv = d["state_variables"];
        state_variables = sv;
    }

    if (d.has("room_types")) {
        Array arr = d["room_types"];
        for (int i = 0; i < arr.size(); i++) {
            Ref<LevelGrammarRoomTypeSpec> rt;
            rt.instantiate();
            rt->from_dictionary(arr[i]);
            room_types.append(rt);
        }
    }

    if (d.has("rules")) {
        Array arr = d["rules"];
        for (int i = 0; i < arr.size(); i++) {
            Ref<LevelGrammarRuleSpec> r;
            r.instantiate();
            r->from_dictionary(arr[i]);
            rules.append(r);
        }
    }
}

// ── JSON Serialization ────────────────────────────────────────────────────

String LevelGrammarSpec::to_json(bool indent) const {
    Dictionary d = to_dictionary();
    String result = JSON::stringify(d, indent ? "    " : "");
    return result;
}

Error LevelGrammarSpec::from_json(const String &json_str,
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

Error LevelGrammarSpec::_from_json_bound(const String &json_str) {
    return from_json(json_str);
}

Error LevelGrammarSpec::save_to_file(const String &path) {
    return ResourceSaver::get_singleton()->save(Ref<Resource>(this), path);
}

Ref<LevelGrammarSpec> LevelGrammarSpec::load_from_json_file(const String &path) {
    Ref<FileAccess> f = FileAccess::open(path, FileAccess::READ);
    if (f.is_null()) {
        UtilityFunctions::printerr("LevelGrammarSpec: Could not open file: ", path);
        return Ref<LevelGrammarSpec>();
    }
    String content = f->get_as_text();
    f->close();

    Ref<LevelGrammarSpec> spec;
    spec.instantiate();
    Error err = spec->from_json(content);
    if (err != OK) {
        UtilityFunctions::printerr("LevelGrammarSpec: JSON parse error in ", path, ": ", err);
        return Ref<LevelGrammarSpec>();
    }
    return spec;
}

} // namespace godot
