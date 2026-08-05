#ifndef LEVEL_GRAMMAR_SPEC_H
#define LEVEL_GRAMMAR_SPEC_H

#include <godot_cpp/classes/resource.hpp>
#include <godot_cpp/variant/array.hpp>
#include <godot_cpp/variant/dictionary.hpp>
#include <godot_cpp/variant/string.hpp>
#include <godot_cpp/variant/typed_array.hpp>
#include <godot_cpp/variant/vector3i.hpp>
#include <godot_cpp/variant/color.hpp>

namespace godot {

// ─────────────────────────────────────────────────────────────────────────────
// LevelGrammarRoomType  —  describes one symbol in the palette
// ─────────────────────────────────────────────────────────────────────────────
class LevelGrammarRoomTypeSpec : public Resource {
    GDCLASS(LevelGrammarRoomTypeSpec, Resource);

private:
    String symbol = "Room";
    Color color = Color(0.4f, 0.6f, 0.9f);
    float weight = 1.0f;
    Vector3i min_size = Vector3i(5, 3, 5);
    Vector3i max_size = Vector3i(10, 6, 10);
    String vox_path;
    bool exclude_from_smoothing = false;
    bool exclude_from_warping = false;

protected:
    static void _bind_methods();

public:
    LevelGrammarRoomTypeSpec();
    virtual ~LevelGrammarRoomTypeSpec() = default;

    // Accessors
    void set_symbol(const String &p_sym);
    String get_symbol() const;
    void set_color(const Color &p_color);
    Color get_color() const;
    void set_weight(float p_w);
    float get_weight() const;
    void set_min_size(const Vector3i &p_size);
    Vector3i get_min_size() const;
    void set_max_size(const Vector3i &p_size);
    Vector3i get_max_size() const;
    void set_vox_path(const String &p_path);
    String get_vox_path() const;
    void set_exclude_from_smoothing(bool p_exclude);
    bool get_exclude_from_smoothing() const;
    void set_exclude_from_warping(bool p_exclude);
    bool get_exclude_from_warping() const;

    // Serialization
    Dictionary to_dictionary() const;
    void from_dictionary(const Dictionary &d);
};


// ─────────────────────────────────────────────────────────────────────────────
// LevelGrammarRuleSpec  —  a single grammar expansion rule
// ─────────────────────────────────────────────────────────────────────────────
class LevelGrammarRuleSpec : public Resource {
    GDCLASS(LevelGrammarRuleSpec, Resource);

private:
    String rule_name = "New Rule";
    String lhs_symbol = "Room";
    float probability = 1.0f;
    String entry_node_id;
    String exit_node_id;

    // Condition
    String condition_var;
    String condition_op = "<";
    float condition_val = 0.0f;

    // Actions: Array of { "var": String, "delta": float }
    Array actions;

    // RHS graph data
    Array rhs_nodes;  // Array of { "id", "symbol", "constraints"?, ... }
    Array rhs_edges;  // Array of { "from", "to", "type" }

protected:
    static void _bind_methods();

public:
    LevelGrammarRuleSpec();
    virtual ~LevelGrammarRuleSpec() = default;

    // Accessors
    void set_rule_name(const String &p_name);
    String get_rule_name() const;
    void set_lhs_symbol(const String &p_sym);
    String get_lhs_symbol() const;
    void set_probability(float p);
    float get_probability() const;
    void set_entry_node_id(const String &p_id);
    String get_entry_node_id() const;
    void set_exit_node_id(const String &p_id);
    String get_exit_node_id() const;
    void set_condition_var(const String &p_var);
    String get_condition_var() const;
    void set_condition_op(const String &p_op);
    String get_condition_op() const;
    void set_condition_val(float p_val);
    float get_condition_val() const;
    void set_actions(const Array &p_actions);
    Array get_actions() const;
    void set_rhs_nodes(const Array &p_nodes);
    Array get_rhs_nodes() const;
    void set_rhs_edges(const Array &p_edges);
    Array get_rhs_edges() const;

    // Serialization
    Dictionary to_dictionary() const;
    void from_dictionary(const Dictionary &d);
};


// ─────────────────────────────────────────────────────────────────────────────
// LevelGrammarSpec  —  top-level grammar definition
//
// This is a C++ Resource that mirrors the GDScript LevelGrammarResource
// property-by-property.  Save it as .tres and the existing
// UnderGenGrammarNode will load and use it transparently.
//
// Builder usage (GDScript / C++):
//
//   var spec = LevelGrammarSpec.new()
//   spec.set_axiom("Start")
//   spec.add_room_type("Entry", Color.CYAN, 1.0, Vector3i(4,3,4), Vector3i(8,5,8))
//   spec.add_rule("Main", "Start", 1.0, [
//       {"id":"a","symbol":"Entry"},
//       {"id":"b","symbol":"Hallway"}
//   ], [
//       {"from":"a","to":"b","type":"corridor"}
//   ])
//   ResourceSaver.save(spec, "res://my_grammar.tres")
//
// JSON usage:
//
//   var spec = LevelGrammarSpec.new()
//   spec.from_json(json_string)
//   ResourceSaver.save(spec, "res://my_grammar.tres")
//
// ─────────────────────────────────────────────────────────────────────────────
class LevelGrammarSpec : public Resource {
    GDCLASS(LevelGrammarSpec, Resource);

private:
    String axiom = "Start";
    Array state_variables;
    TypedArray<LevelGrammarRoomTypeSpec> room_types;
    TypedArray<LevelGrammarRuleSpec> rules;

protected:
    static void _bind_methods();

public:
    LevelGrammarSpec();
    virtual ~LevelGrammarSpec() = default;

    // Accessors
    void set_axiom(const String &p_axiom);
    String get_axiom() const;
    void set_state_variables(const Array &p_vars);
    Array get_state_variables() const;
    void set_room_types(const TypedArray<LevelGrammarRoomTypeSpec> &p_types);
    TypedArray<LevelGrammarRoomTypeSpec> get_room_types() const;
    void set_rules(const TypedArray<LevelGrammarRuleSpec> &p_rules);
    TypedArray<LevelGrammarRuleSpec> get_rules() const;

    // ── Builder API ───────────────────────────────────────────────────────
    /// Add a room type with one call. Returns the created spec for chaining.
    Ref<LevelGrammarRoomTypeSpec> add_room_type(
        const String &symbol,
        const Color &color = Color(0.4f, 0.6f, 0.9f),
        float weight = 1.0f,
        const Vector3i &min_size = Vector3i(5, 3, 5),
        const Vector3i &max_size = Vector3i(10, 6, 10),
        const String &vox_path = "");

    /// Convenience: add room type with size defaults, optionally a vox.
    Ref<LevelGrammarRoomTypeSpec> add_room_type_vox(
        const String &symbol, const String &vox_path,
        const Color &color = Color(0.4f, 0.6f, 0.9f));

    /// Add a grammar rule with RHS nodes/edges as arrays of dictionaries.
    /// Returns the created rule spec for chaining.
    Ref<LevelGrammarRuleSpec> add_rule(
        const String &rule_name,
        const String &lhs_symbol,
        float probability,
        const Array &rhs_nodes,
        const Array &rhs_edges,
        const String &entry_node_id = "",
        const String &exit_node_id = "");

    /// Add a rule with a condition.
    Ref<LevelGrammarRuleSpec> add_rule_conditional(
        const String &rule_name,
        const String &lhs_symbol,
        float probability,
        const String &condition_var,
        const String &condition_op,
        float condition_val,
        const Array &actions,
        const Array &rhs_nodes,
        const Array &rhs_edges,
        const String &entry_node_id = "",
        const String &exit_node_id = "");

    // ── Serialization ─────────────────────────────────────────────────────
    /// Serialize the complete grammar to a Godot Dictionary.
    Dictionary to_dictionary() const;

    /// Populate the grammar from a Godot Dictionary (same structure).
    void from_dictionary(const Dictionary &d);

    /// Serialize to a JSON string.
    String to_json(bool indent = true) const;

    /// Populate from a JSON string. Returns OK/FAIL via Error enum.
    /// On failure, err_line / err_message are filled.
    Error from_json(const String &json_str,
                    int *r_err_line = nullptr,
                    String *r_err_message = nullptr);

    /// ClassDB-bindable overload (no pointer params).
    Error _from_json_bound(const String &json_str);

    /// Save this resource to disk as a .tres that the GDScript editor can open.
    Error save_to_file(const String &path);

    /// Static factory: create and populate from a JSON file path.
    static Ref<LevelGrammarSpec> load_from_json_file(const String &path);
};

} // namespace godot

#endif // LEVEL_GRAMMAR_SPEC_H
