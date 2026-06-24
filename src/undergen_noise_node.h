#ifndef UNDERGEN_NOISE_NODE_H
#define UNDERGEN_NOISE_NODE_H

#include "undergen_node.h"
#include <godot_cpp/classes/fast_noise_lite.hpp>

namespace godot {

class UnderGenNoiseNode : public UnderGenNode {
    GDCLASS(UnderGenNoiseNode, UnderGenNode);

private:
    float noise_scale = 50.0f;
    float noise_intensity = 0.5f;
    float noise_frequency = 0.02f;
    int noise_seed = 1337;
    Ref<FastNoiseLite> noise_generator;

protected:
    static void _bind_methods();

public:
    UnderGenNoiseNode();
    virtual ~UnderGenNoiseNode();

    void set_noise_scale(float p_scale);
    float get_noise_scale() const;
    void set_noise_intensity(float p_intensity);
    float get_noise_intensity() const;
    void set_noise_frequency(float p_freq);
    float get_noise_frequency() const;
    void set_noise_seed(int p_seed);
    int get_noise_seed() const;
    void set_noise_generator(const Ref<FastNoiseLite> &p_noise);
    Ref<FastNoiseLite> get_noise_generator() const;

    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
};

} // namespace godot

#endif // UNDERGEN_NOISE_NODE_H
