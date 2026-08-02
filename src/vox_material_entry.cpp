#include "vox_material_entry.h"

#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/core/math.hpp>

namespace godot {

VoxMaterialEntry::VoxMaterialEntry() {
    thickness = 1.0f;
    stepped = true;
}
VoxMaterialEntry::~VoxMaterialEntry() {}

void VoxMaterialEntry::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_palette_index", "index"), &VoxMaterialEntry::set_palette_index);
    ClassDB::bind_method(D_METHOD("get_palette_index"), &VoxMaterialEntry::get_palette_index);
    ClassDB::bind_method(D_METHOD("set_material_id", "id"), &VoxMaterialEntry::set_material_id);
    ClassDB::bind_method(D_METHOD("get_material_id"), &VoxMaterialEntry::get_material_id);
    ClassDB::bind_method(D_METHOD("set_thickness", "thickness"), &VoxMaterialEntry::set_thickness);
    ClassDB::bind_method(D_METHOD("get_thickness"), &VoxMaterialEntry::get_thickness);
    ClassDB::bind_method(D_METHOD("set_stepped", "stepped"), &VoxMaterialEntry::set_stepped);
    ClassDB::bind_method(D_METHOD("is_stepped"), &VoxMaterialEntry::is_stepped);

    ADD_PROPERTY(PropertyInfo(Variant::INT, "palette_index", PROPERTY_HINT_RANGE, "0,255,1"), "set_palette_index", "get_palette_index");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "material_id", PROPERTY_HINT_RANGE, "0,255,1"), "set_material_id", "get_material_id");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "thickness", PROPERTY_HINT_RANGE, "0.05,5.0,0.05"), "set_thickness", "get_thickness");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "stepped"), "set_stepped", "is_stepped");
}

void VoxMaterialEntry::set_palette_index(int p_index) {
    palette_index = Math::clamp(p_index, 0, 255);
    emit_changed();
}

int VoxMaterialEntry::get_palette_index() const {
    return palette_index;
}

void VoxMaterialEntry::set_material_id(int p_id) {
    material_id = Math::clamp(p_id, 0, 255);
    emit_changed();
}

int VoxMaterialEntry::get_material_id() const {
    return material_id;
}

void VoxMaterialEntry::set_thickness(float p_thickness) {
    thickness = Math::clamp(p_thickness, 0.05f, 5.0f);
    emit_changed();
}

float VoxMaterialEntry::get_thickness() const {
    return thickness;
}

void VoxMaterialEntry::set_stepped(bool p_stepped) {
    stepped = p_stepped;
    emit_changed();
}

bool VoxMaterialEntry::is_stepped() const {
    return stepped;
}

} // namespace godot
