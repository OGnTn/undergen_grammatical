#ifndef UNDERGEN_POINT_SET_H
#define UNDERGEN_POINT_SET_H

#include <godot_cpp/classes/ref_counted.hpp>
#include <godot_cpp/variant/transform3d.hpp>
#include <godot_cpp/variant/aabb.hpp>
#include <godot_cpp/variant/color.hpp>
#include <godot_cpp/variant/dictionary.hpp>
#include <vector>

namespace godot {

class UnderGenPointSet : public RefCounted {
    GDCLASS(UnderGenPointSet, RefCounted);

public:
    struct PCGPoint {
        Transform3D transform;
        float density = 1.0f;
        AABB bounds;
        Color color;
        Dictionary attributes;
    };

private:
    std::vector<PCGPoint> points;

protected:
    static void _bind_methods();

public:
    UnderGenPointSet();
    virtual ~UnderGenPointSet();

    void add_point(const Transform3D &transform, float density, const AABB &bounds, const Color &color, const Dictionary &attributes);
    void add_raw_point(const Transform3D &transform, float density, const Dictionary &attributes); // Simpler helper
    void remove_point(int idx);
    int get_point_count() const;

    Transform3D get_point_transform(int idx) const;
    float get_point_density(int idx) const;
    AABB get_point_bounds(int idx) const;
    Color get_point_color(int idx) const;
    Dictionary get_point_attributes(int idx) const;

    void clear();
    Array get_points() const;

    // Direct access for C++ nodes
    std::vector<PCGPoint>& get_points_raw() { return points; }
    const std::vector<PCGPoint>& get_points_raw() const { return points; }
};

} // namespace godot

#endif // UNDERGEN_POINT_SET_H
