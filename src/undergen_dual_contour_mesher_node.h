// src/undergen_dual_contour_mesher_node.h
#ifndef UNDERGEN_DUAL_CONTOUR_MESHER_NODE_H
#define UNDERGEN_DUAL_CONTOUR_MESHER_NODE_H

#include "undergen_mesher_node.h"

namespace godot {

class UnderGenDualContourMesherNode : public UnderGenMesherNode {
    GDCLASS(UnderGenDualContourMesherNode, UnderGenMesherNode);

private:
    bool use_qef = true;
    float qef_regularization = 1e-4f;

protected:
    static void _bind_methods();

public:
    UnderGenDualContourMesherNode();
    virtual ~UnderGenDualContourMesherNode();

    void set_use_qef(bool p_use);
    bool get_use_qef() const;

    void set_qef_regularization(float p_reg);
    float get_qef_regularization() const;

    virtual void execute_with_parent(const Dictionary &inputs, Dictionary &outputs, Node3D* parent_node) override;
};

} // namespace godot

#endif // UNDERGEN_DUAL_CONTOUR_MESHER_NODE_H
