#ifndef UNDERGEN_SPLINE_CARVER_NODE_H
#define UNDERGEN_SPLINE_CARVER_NODE_H

#include "undergen_node.h"
#include <godot_cpp/classes/curve3d.hpp>
#include <godot_cpp/variant/packed_vector2_array.hpp>
#include <godot_cpp/variant/packed_vector3_array.hpp>
#include <godot_cpp/variant/string.hpp>
#include <godot_cpp/variant/vector3i.hpp>

namespace godot {

class DensityGrid;

class UnderGenSplineCarverNode : public UnderGenNode {
    GDCLASS(UnderGenSplineCarverNode, UnderGenNode);

public:
    enum ProfileType {
        PROFILE_CIRCLE = 0,
        PROFILE_HORSESHOE = 1,
        PROFILE_RECTANGLE = 2,
        PROFILE_GOTHIC_ARCH = 3,
        PROFILE_CUSTOM = 4
    };

    enum CarveMode {
        CARVE_SUBTRACT = 0,
        CARVE_ADD = 1
    };

private:
    ProfileType profile_type = PROFILE_HORSESHOE;
    CarveMode carve_mode = CARVE_SUBTRACT;

    float width = 6.0f;
    float height = 5.0f;
    float wall_height = 3.0f;
    float sample_step = 0.5f;
    float profile_rotation_deg = 0.0f;

    PackedVector2Array custom_profile;
    String curve_resource_path = "";
    Ref<Curve3D> curve;

    bool stamp_materials = false;
    int floor_material_id = 0;
    int wall_material_id = 0;
    int ceiling_material_id = 0;

    float _evaluate_2d_sdf(float u, float v) const;
    void _carve_single_curve(Ref<DensityGrid> grid, const Ref<Curve3D> &active_curve);

protected:
    static void _bind_methods();

public:
    UnderGenSplineCarverNode();
    virtual ~UnderGenSplineCarverNode();

    void set_profile_type(ProfileType p_type);
    ProfileType get_profile_type() const;

    void set_carve_mode(CarveMode p_mode);
    CarveMode get_carve_mode() const;

    void set_width(float p_width);
    float get_width() const;

    void set_height(float p_height);
    float get_height() const;

    void set_wall_height(float p_wheight);
    float get_wall_height() const;

    void set_sample_step(float p_step);
    float get_sample_step() const;

    void set_profile_rotation_deg(float p_deg);
    float get_profile_rotation_deg() const;

    void set_custom_profile(const PackedVector2Array &p_profile);
    PackedVector2Array get_custom_profile() const;

    void set_curve_resource_path(const String &p_path);
    String get_curve_resource_path() const;

    void set_curve(const Ref<Curve3D> &p_curve);
    Ref<Curve3D> get_curve() const;

    void set_stamp_materials(bool p_stamp);
    bool get_stamp_materials() const;

    void set_floor_material_id(int p_id);
    int get_floor_material_id() const;

    void set_wall_material_id(int p_id);
    int get_wall_material_id() const;

    void set_ceiling_material_id(int p_id);
    int get_ceiling_material_id() const;

    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
};

} // namespace godot

VARIANT_ENUM_CAST(godot::UnderGenSplineCarverNode::ProfileType);
VARIANT_ENUM_CAST(godot::UnderGenSplineCarverNode::CarveMode);

#endif // UNDERGEN_SPLINE_CARVER_NODE_H
