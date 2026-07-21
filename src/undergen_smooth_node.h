#ifndef UNDERGEN_SMOOTH_NODE_H
#define UNDERGEN_SMOOTH_NODE_H

#include "undergen_node.h"

namespace godot {

class UnderGenSmoothNode : public UnderGenNode {
    GDCLASS(UnderGenSmoothNode, UnderGenNode);

private:
    int smoothing_strength = 1;
    bool enforce_solid_boundaries = true;
    static constexpr float WORLD_SOLID_VALUE = 1.0f;

protected:
    static void _bind_methods();

public:
    UnderGenSmoothNode();
    virtual ~UnderGenSmoothNode();

    void set_smoothing_strength(int p_strength);
    int get_smoothing_strength() const;
    void set_enforce_solid_boundaries(bool p_enforce);
    bool get_enforce_solid_boundaries() const;

    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
};

} // namespace godot

#endif // UNDERGEN_SMOOTH_NODE_H
