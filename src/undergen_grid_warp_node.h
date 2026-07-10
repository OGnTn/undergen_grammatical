#ifndef UNDERGEN_GRID_WARP_NODE_H
#define UNDERGEN_GRID_WARP_NODE_H

#include "undergen_node.h"
#include <godot_cpp/classes/fast_noise_lite.hpp>

namespace godot {

class UnderGenGridWarpNode : public UnderGenNode {
    GDCLASS(UnderGenGridWarpNode, UnderGenNode);

private:
    float warp_amplitude = 4.0f;
    float noise_frequency = 0.03f;
    int noise_seed = 1337;
    Ref<FastNoiseLite> noise_generator;

protected:
    static void _bind_methods();

public:
    UnderGenGridWarpNode();
    virtual ~UnderGenGridWarpNode();

    void set_warp_amplitude(float p_amplitude);
    float get_warp_amplitude() const;
    void set_noise_frequency(float p_freq);
    float get_noise_frequency() const;
    void set_noise_seed(int p_seed);
    int get_noise_seed() const;
    void set_noise_generator(const Ref<FastNoiseLite> &p_noise);
    Ref<FastNoiseLite> get_noise_generator() const;

    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
};

} // namespace godot

#endif // UNDERGEN_GRID_WARP_NODE_H
