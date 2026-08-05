#include "undergen_node.h"

namespace godot {

UnderGenNode::UnderGenNode() : editor_position(Vector2(0, 0)) {}
UnderGenNode::~UnderGenNode() {}

void UnderGenNode::_bind_methods() {
    ClassDB::bind_method(D_METHOD("execute", "inputs"), &UnderGenNode::execute_to_dictionary);
    ClassDB::bind_method(D_METHOD("set_editor_position", "editor_position"), &UnderGenNode::set_editor_position);
    ClassDB::bind_method(D_METHOD("get_editor_position"), &UnderGenNode::get_editor_position);

    ADD_PROPERTY(PropertyInfo(Variant::VECTOR2, "editor_position"), "set_editor_position", "get_editor_position");
}

void UnderGenNode::set_editor_position(const Vector2 &p_pos) {
    editor_position = p_pos;
}

Vector2 UnderGenNode::get_editor_position() const {
    return editor_position;
}

void UnderGenNode::execute(const Dictionary &inputs, Dictionary &outputs) {
    _execute(inputs, outputs);
}

Dictionary UnderGenNode::execute_to_dictionary(const Dictionary &inputs) {
    Dictionary outputs;
    execute(inputs, outputs);
    return outputs;
}

void UnderGenNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    // Default implementation does nothing.
}

Dictionary UnderGenNode::get_pipeline_input_defaults(const Dictionary &global_inputs) const {
    return Dictionary();
}

} // namespace godot
