#include "undergen_point_set.h"

namespace godot {

UnderGenPointSet::UnderGenPointSet() {}
UnderGenPointSet::~UnderGenPointSet() {}

void UnderGenPointSet::_bind_methods() {
    ClassDB::bind_method(D_METHOD("add_point", "transform", "density", "bounds", "color", "attributes"), &UnderGenPointSet::add_point);
    ClassDB::bind_method(D_METHOD("add_raw_point", "transform", "density", "attributes"), &UnderGenPointSet::add_raw_point);
    ClassDB::bind_method(D_METHOD("remove_point", "idx"), &UnderGenPointSet::remove_point);
    ClassDB::bind_method(D_METHOD("get_point_count"), &UnderGenPointSet::get_point_count);

    ClassDB::bind_method(D_METHOD("get_point_transform", "idx"), &UnderGenPointSet::get_point_transform);
    ClassDB::bind_method(D_METHOD("get_point_density", "idx"), &UnderGenPointSet::get_point_density);
    ClassDB::bind_method(D_METHOD("get_point_bounds", "idx"), &UnderGenPointSet::get_point_bounds);
    ClassDB::bind_method(D_METHOD("get_point_color", "idx"), &UnderGenPointSet::get_point_color);
    ClassDB::bind_method(D_METHOD("get_point_attributes", "idx"), &UnderGenPointSet::get_point_attributes);

    ClassDB::bind_method(D_METHOD("clear"), &UnderGenPointSet::clear);
    ClassDB::bind_method(D_METHOD("get_points"), &UnderGenPointSet::get_points);
}

void UnderGenPointSet::add_point(const Transform3D &transform, float density, const AABB &bounds, const Color &color, const Dictionary &attributes) {
    PCGPoint p;
    p.transform = transform;
    p.density = density;
    p.bounds = bounds;
    p.color = color;
    p.attributes = attributes.duplicate();
    points.push_back(p);
}

void UnderGenPointSet::add_raw_point(const Transform3D &transform, float density, const Dictionary &attributes) {
    PCGPoint p;
    p.transform = transform;
    p.density = density;
    p.bounds = AABB(Vector3(-0.5, -0.5, -0.5), Vector3(1.0, 1.0, 1.0));
    p.color = Color(1, 1, 1, 1);
    p.attributes = attributes.duplicate();
    points.push_back(p);
}

void UnderGenPointSet::remove_point(int idx) {
    if (idx >= 0 && idx < (int)points.size()) {
        points.erase(points.begin() + idx);
    }
}

int UnderGenPointSet::get_point_count() const {
    return (int)points.size();
}

Transform3D UnderGenPointSet::get_point_transform(int idx) const {
    if (idx >= 0 && idx < (int)points.size()) {
        return points[idx].transform;
    }
    return Transform3D();
}

float UnderGenPointSet::get_point_density(int idx) const {
    if (idx >= 0 && idx < (int)points.size()) {
        return points[idx].density;
    }
    return 0.0f;
}

AABB UnderGenPointSet::get_point_bounds(int idx) const {
    if (idx >= 0 && idx < (int)points.size()) {
        return points[idx].bounds;
    }
    return AABB();
}

Color UnderGenPointSet::get_point_color(int idx) const {
    if (idx >= 0 && idx < (int)points.size()) {
        return points[idx].color;
    }
    return Color();
}

Dictionary UnderGenPointSet::get_point_attributes(int idx) const {
    if (idx >= 0 && idx < (int)points.size()) {
        return points[idx].attributes;
    }
    return Dictionary();
}

void UnderGenPointSet::clear() {
    points.clear();
}

Array UnderGenPointSet::get_points() const {
    Array arr;
    for (const auto &p : points) {
        Dictionary d;
        d["transform"] = p.transform;
        d["density"] = p.density;
        d["bounds"] = p.bounds;
        d["color"] = p.color;
        d["attributes"] = p.attributes;
        arr.append(d);
    }
    return arr;
}

} // namespace godot
