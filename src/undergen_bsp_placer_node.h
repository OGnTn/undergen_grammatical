#ifndef UNDERGEN_BSP_PLACER_NODE_H
#define UNDERGEN_BSP_PLACER_NODE_H

#include "undergen_node.h"
#include "density_grid.h"
#include "room_generator.h"
#include <godot_cpp/classes/random_number_generator.hpp>

namespace godot {

class UnderGenBSPPlacerNode : public UnderGenNode {
    GDCLASS(UnderGenBSPPlacerNode, UnderGenNode);

private:
    int grid_size_x = 48;
    int grid_size_y = 16;
    int grid_size_z = 48;
    Ref<RandomNumberGenerator> rng;

protected:
    static void _bind_methods();

public:
    UnderGenBSPPlacerNode();
    virtual ~UnderGenBSPPlacerNode();

    void set_grid_size_x(int p_x);
    int get_grid_size_x() const;
    void set_grid_size_y(int p_y);
    int get_grid_size_y() const;
    void set_grid_size_z(int p_z);
    int get_grid_size_z() const;

    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
};

} // namespace godot

#endif // UNDERGEN_BSP_PLACER_NODE_H
