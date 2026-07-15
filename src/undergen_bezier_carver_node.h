#ifndef UNDERGEN_BEZIER_CARVER_NODE_H
#define UNDERGEN_BEZIER_CARVER_NODE_H

#include "undergen_node.h"
#include "path_carver.h"
#include <godot_cpp/classes/fast_noise_lite.hpp>
#include <godot_cpp/classes/random_number_generator.hpp>

namespace godot {

class UnderGenBezierCarverNode : public UnderGenNode {
    GDCLASS(UnderGenBezierCarverNode, UnderGenNode);

private:
    // Brush
    int path_brush_min_radius = 2;
    int path_brush_max_radius = 4;
    bool use_square_brush = false;

    // Bezier
    int path_segments = 1;
    float path_bend_factor = 0.4f;
    float path_wobble_magnitude = 0.0f;
    float path_wobble_frequency = 0.2f;
    bool connect_from_ground_level = false;
    bool dont_carve_inside_rooms = false;

    // Organic Cave Shape
    float cave_ruggedness = 1.0f;
    float cave_floor_ruggedness = 0.0f;
    float cave_ceiling_ruggedness = 0.0f;
    float cave_width_noise = 0.0f;
    float floor_flattening = 0.0f;
    float overhang_openness = 0.0f;

    Ref<RandomNumberGenerator> rng;
    Ref<FastNoiseLite> wobble_noise;

protected:
    static void _bind_methods();

public:
    UnderGenBezierCarverNode();
    virtual ~UnderGenBezierCarverNode();

    // Brush
    void set_path_brush_min_radius(int p_radius);
    int get_path_brush_min_radius() const;
    void set_path_brush_max_radius(int p_radius);
    int get_path_brush_max_radius() const;
    void set_use_square_brush(bool p_enabled);
    bool get_use_square_brush() const;

    // Bezier
    void set_path_segments(int p_segments);
    int get_path_segments() const;
    void set_path_bend_factor(float p_factor);
    float get_path_bend_factor() const;
    void set_path_wobble_magnitude(float p_magnitude);
    float get_path_wobble_magnitude() const;
    void set_path_wobble_frequency(float p_frequency);
    float get_path_wobble_frequency() const;
    void set_connect_from_ground_level(bool p_enabled);
    bool get_connect_from_ground_level() const;
    void set_dont_carve_inside_rooms(bool p_enabled);
    bool get_dont_carve_inside_rooms() const;

    // Organic Cave Shape
    void set_cave_ruggedness(float p_value);
    float get_cave_ruggedness() const;
    void set_cave_floor_ruggedness(float p_value);
    float get_cave_floor_ruggedness() const;
    void set_cave_ceiling_ruggedness(float p_value);
    float get_cave_ceiling_ruggedness() const;
    void set_cave_width_noise(float p_value);
    float get_cave_width_noise() const;
    void set_floor_flattening(float p_value);
    float get_floor_flattening() const;
    void set_overhang_openness(float p_value);
    float get_overhang_openness() const;

    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
};

} // namespace godot

#endif // UNDERGEN_BEZIER_CARVER_NODE_H
