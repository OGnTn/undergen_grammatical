#ifndef UNDERGEN_OUTDOOR_PLACER_NODE_H
#define UNDERGEN_OUTDOOR_PLACER_NODE_H

#include "undergen_node.h"
#include "density_grid.h"
#include <godot_cpp/classes/random_number_generator.hpp>
#include <godot_cpp/classes/fast_noise_lite.hpp>

namespace godot {

class UnderGenOutdoorPlacerNode : public UnderGenNode {
    GDCLASS(UnderGenOutdoorPlacerNode, UnderGenNode);

private:
    int grid_size_x = 48;
    int grid_size_y = 16;
    int grid_size_z = 48;
    int margin_x = 1;
    int margin_y = 1;
    int margin_z = 1;
    float surface_threshold = 0.0f;
    float spread_ratio = 0.5f;
    int base_height = 4;
    int mountain_height = 12;
    float slope_width = 4.0f;
    float path_width = 3.0f;
    float noise_intensity = 3.0f;
    Ref<FastNoiseLite> noise_generator;
    Ref<RandomNumberGenerator> rng;

protected:
    static void _bind_methods();

public:
    UnderGenOutdoorPlacerNode();
    virtual ~UnderGenOutdoorPlacerNode();

    void set_grid_size_x(int p_x);
    int get_grid_size_x() const;
    void set_grid_size_y(int p_y);
    int get_grid_size_y() const;
    void set_grid_size_z(int p_z);
    int get_grid_size_z() const;
    void set_surface_threshold(float p_threshold);
    float get_surface_threshold() const;

    void set_margin_x(int p_x);
    int get_margin_x() const;
    void set_margin_y(int p_y);
    int get_margin_y() const;
    void set_margin_z(int p_z);
    int get_margin_z() const;

    void set_spread_ratio(float p_ratio);
    float get_spread_ratio() const;

    void set_base_height(int p_height);
    int get_base_height() const;
    void set_mountain_height(int p_height);
    int get_mountain_height() const;
    void set_slope_width(float p_width);
    float get_slope_width() const;
    void set_path_width(float p_width);
    float get_path_width() const;
    void set_noise_intensity(float p_intensity);
    float get_noise_intensity() const;
    void set_noise_generator(const Ref<FastNoiseLite> &p_noise);
    Ref<FastNoiseLite> get_noise_generator() const;

    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
};

} // namespace godot

#endif // UNDERGEN_OUTDOOR_PLACER_NODE_H
