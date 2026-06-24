#ifndef UNDERGEN_NODE_H
#define UNDERGEN_NODE_H

#include <godot_cpp/classes/resource.hpp>
#include <godot_cpp/variant/dictionary.hpp>
#include <godot_cpp/variant/vector2.hpp>

namespace godot {

class UnderGenNode : public Resource {
    GDCLASS(UnderGenNode, Resource);

private:
    Vector2 editor_position;

protected:
    static void _bind_methods();

public:
    UnderGenNode();
    virtual ~UnderGenNode();

    void set_editor_position(const Vector2 &p_pos);
    Vector2 get_editor_position() const;

    virtual void execute(const Dictionary &inputs, Dictionary &outputs);
    virtual void _execute(const Dictionary &inputs, Dictionary &outputs);
};

} // namespace godot

#endif // UNDERGEN_NODE_H
