// src/undergen_grammar_node.h
#ifndef UNDERGEN_GRAMMAR_NODE_H
#define UNDERGEN_GRAMMAR_NODE_H

#include "undergen_node.h"
#include <godot_cpp/classes/random_number_generator.hpp>
#include <vector>
#include <map>
#include <string>

namespace godot {

/// UnderGenGrammarNode
/// ---
/// Pipeline node that runs grammar expansion on a LevelGrammarResource.
/// Input  port 0: Seed (int)
/// Output port 0: Logical Graph (Dictionary { "nodes": Array, "edges": Array })
///
/// Connect its output directly to UnderGenBSPPlacerNode input port 1.
class UnderGenGrammarNode : public UnderGenNode {
    GDCLASS(UnderGenGrammarNode, UnderGenNode);

private:
    String grammar_resource_path;
    int iterations = 4;
    int max_nodes  = 100;

    Ref<RandomNumberGenerator> rng;

    // --------------- internal helpers ---------------
    bool _eval_condition(
        const String& var, const String& op, double val,
        const Dictionary& state) const;

    // Returns a rule resource chosen by weighted probability,
    // or an invalid Ref if no matching rule passes its condition.
    Ref<Resource> _pick_rule(const Array& candidates, const Dictionary& state);

    // Enrich a raw node dictionary with room-type size/vox data.
    Dictionary _make_node(
        const String& id, const String& symbol,
        const Dictionary& constraints,
        const std::map<String, Dictionary>& rt_map) const;

protected:
    static void _bind_methods();

public:
    UnderGenGrammarNode();
    virtual ~UnderGenGrammarNode() = default;

    void set_grammar_resource_path(const String& p_path);
    String get_grammar_resource_path() const;
    void set_iterations(int p_iters);
    int  get_iterations() const;
    void set_max_nodes(int p_max);
    int  get_max_nodes() const;

    virtual void _execute(const Dictionary& inputs, Dictionary& outputs) override;
    virtual Dictionary get_pipeline_input_defaults(const Dictionary &global_inputs) const override;
};

} // namespace godot
#endif // UNDERGEN_GRAMMAR_NODE_H
