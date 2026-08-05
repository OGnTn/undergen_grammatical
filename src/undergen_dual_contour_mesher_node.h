// src/undergen_dual_contour_mesher_node.h
#ifndef UNDERGEN_DUAL_CONTOUR_MESHER_NODE_H
#define UNDERGEN_DUAL_CONTOUR_MESHER_NODE_H

#include "undergen_mesher_node.h"

namespace godot {

class UnderGenDualContourMesherNode : public UnderGenMesherNode {
    GDCLASS(UnderGenDualContourMesherNode, UnderGenMesherNode);

private:
    bool use_qef = true;
    float qef_regularization = 0.01f;
    bool stepped_transitions = true;
    bool adaptive_mesh = false;
    float curvature_threshold = 0.05f;

protected:
    static void _bind_methods();

public:
    UnderGenDualContourMesherNode();
    virtual ~UnderGenDualContourMesherNode();

    void set_use_qef(bool p_use);
    bool get_use_qef() const;

    void set_qef_regularization(float p_reg);
    float get_qef_regularization() const;

    void set_stepped_transitions(bool p_stepped);
    bool get_stepped_transitions() const;

    void set_adaptive_mesh(bool p_adaptive);
    bool get_adaptive_mesh() const;

    void set_curvature_threshold(float p_threshold);
    float get_curvature_threshold() const;

    virtual void execute_with_parent(const Dictionary &inputs, Dictionary &outputs, Node3D* parent_node) override;
};

} // namespace godot

#endif // UNDERGEN_DUAL_CONTOUR_MESHER_NODE_H
