#include "undergen_noise_node.h"
#include "density_grid.h"
#include <godot_cpp/variant/utility_functions.hpp>

namespace godot {

UnderGenNoiseNode::UnderGenNoiseNode() {
    noise_generator.instantiate();
}

UnderGenNoiseNode::~UnderGenNoiseNode() {}

void UnderGenNoiseNode::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_noise_scale", "scale"), &UnderGenNoiseNode::set_noise_scale);
    ClassDB::bind_method(D_METHOD("get_noise_scale"), &UnderGenNoiseNode::get_noise_scale);
    ClassDB::bind_method(D_METHOD("set_noise_intensity", "intensity"), &UnderGenNoiseNode::set_noise_intensity);
    ClassDB::bind_method(D_METHOD("get_noise_intensity"), &UnderGenNoiseNode::get_noise_intensity);
    ClassDB::bind_method(D_METHOD("set_noise_frequency", "frequency"), &UnderGenNoiseNode::set_noise_frequency);
    ClassDB::bind_method(D_METHOD("get_noise_frequency"), &UnderGenNoiseNode::get_noise_frequency);
    ClassDB::bind_method(D_METHOD("set_noise_seed", "seed"), &UnderGenNoiseNode::set_noise_seed);
    ClassDB::bind_method(D_METHOD("get_noise_seed"), &UnderGenNoiseNode::get_noise_seed);
    ClassDB::bind_method(D_METHOD("set_noise_generator", "noise"), &UnderGenNoiseNode::set_noise_generator);
    ClassDB::bind_method(D_METHOD("get_noise_generator"), &UnderGenNoiseNode::get_noise_generator);

    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "noise_scale"), "set_noise_scale", "get_noise_scale");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "noise_intensity"), "set_noise_intensity", "get_noise_intensity");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "noise_frequency"), "set_noise_frequency", "get_noise_frequency");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "noise_seed"), "set_noise_seed", "get_noise_seed");
    ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "noise_generator", PROPERTY_HINT_RESOURCE_TYPE, "FastNoiseLite"), "set_noise_generator", "get_noise_generator");
}

void UnderGenNoiseNode::set_noise_scale(float p_scale) { noise_scale = p_scale; }
float UnderGenNoiseNode::get_noise_scale() const { return noise_scale; }
void UnderGenNoiseNode::set_noise_intensity(float p_intensity) { noise_intensity = p_intensity; }
float UnderGenNoiseNode::get_noise_intensity() const { return noise_intensity; }
void UnderGenNoiseNode::set_noise_frequency(float p_freq) { noise_frequency = p_freq; }
float UnderGenNoiseNode::get_noise_frequency() const { return noise_frequency; }
void UnderGenNoiseNode::set_noise_seed(int p_seed) { noise_seed = p_seed; }
int UnderGenNoiseNode::get_noise_seed() const { return noise_seed; }
void UnderGenNoiseNode::set_noise_generator(const Ref<FastNoiseLite> &p_noise) { noise_generator = p_noise; }
Ref<FastNoiseLite> UnderGenNoiseNode::get_noise_generator() const { return noise_generator; }

void UnderGenNoiseNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    Dictionary context = inputs.get(0, Dictionary());
    if (context.is_empty()) return;

    Ref<DensityGrid> grid = context.get("grid", Ref<DensityGrid>());
    if (grid.is_null() || !noise_generator.is_valid()) return;

    noise_generator->set_seed(noise_seed);
    noise_generator->set_frequency(noise_frequency);

    float inv_noise_scale = (noise_scale > 0.0001f) ? 1.0f / noise_scale : 1.0f;
    int gsx = grid->get_grid_size_x();
    int gsy = grid->get_grid_size_y();
    int gsz = grid->get_grid_size_z();
    float surf_thresh = grid->get_surface_threshold();

    for (int z = 0; z < gsz; ++z) {
        for (int y = 0; y < gsy; ++y) {
            for (int x = 0; x < gsx; ++x) {
                Vector3i pos(x, y, z);
                float noise_val = noise_generator->get_noise_3d(x * inv_noise_scale, y * inv_noise_scale, z * inv_noise_scale);
                float current_cell_value = grid->get_cell(pos, 1.0f);
                float new_density;
                if (current_cell_value < surf_thresh) {
                    float carving_noise = Math::min(noise_val, 0.0f);
                    new_density = current_cell_value + carving_noise * noise_intensity;
                } else {
                    new_density = current_cell_value + noise_val * noise_intensity;
                }
                grid->set_cell(pos, new_density);
            }
        }
    }

    outputs[0] = context;
}

} // namespace godot
