#ifndef UNDERGEN_ASTAR_CARVER_NODE_H
#define UNDERGEN_ASTAR_CARVER_NODE_H

#include "undergen_node.h"
#include "path_carver.h"
#include <godot_cpp/classes/fast_noise_lite.hpp>
#include <godot_cpp/classes/random_number_generator.hpp>

namespace godot {

class UnderGenAStarCarverNode : public UnderGenNode {
    GDCLASS(UnderGenAStarCarverNode, UnderGenNode);

private:
    int path_brush_min_radius = 2;
    int path_brush_max_radius = 4;
    bool use_square_brush = false;
    float vertical_movement_cost_multiplier = 2.0f;
    bool connect_from_ground_level = false;

    Ref<RandomNumberGenerator> rng;
    Ref<FastNoiseLite> wobble_noise;

protected:
    static void _bind_methods();

public:
    UnderGenAStarCarverNode();
    virtual ~UnderGenAStarCarverNode();

    void set_path_brush_min_radius(int p_radius);
    int get_path_brush_min_radius() const;
    void set_path_brush_max_radius(int p_radius);
    int get_path_brush_max_radius() const;
    void set_use_square_brush(bool p_enabled);
    bool get_use_square_brush() const;
    void set_vertical_movement_cost_multiplier(float p_mult);
    float get_vertical_movement_cost_multiplier() const;
    void set_connect_from_ground_level(bool p_enabled);
    bool get_connect_from_ground_level() const;

    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
};

} // namespace godot

#endif // UNDERGEN_ASTAR_CARVER_NODE_H
