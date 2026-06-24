#include "zone_material_entry.h"

#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/core/math.hpp>

namespace godot {

ZoneMaterialEntry::ZoneMaterialEntry() {}
ZoneMaterialEntry::~ZoneMaterialEntry() {}

void ZoneMaterialEntry::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_zone_name", "name"), &ZoneMaterialEntry::set_zone_name);
    ClassDB::bind_method(D_METHOD("get_zone_name"), &ZoneMaterialEntry::get_zone_name);
    ClassDB::bind_method(D_METHOD("set_material_id", "id"), &ZoneMaterialEntry::set_material_id);
    ClassDB::bind_method(D_METHOD("get_material_id"), &ZoneMaterialEntry::get_material_id);

    ADD_PROPERTY(PropertyInfo(Variant::STRING, "zone_name"), "set_zone_name", "get_zone_name");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "material_id", PROPERTY_HINT_RANGE, "0,255,1"), "set_material_id", "get_material_id");
}

void ZoneMaterialEntry::set_zone_name(const String &p_name) { zone_name = p_name; }
String ZoneMaterialEntry::get_zone_name() const { return zone_name; }
void ZoneMaterialEntry::set_material_id(int p_id) { material_id = Math::clamp(p_id, 0, 255); }
int ZoneMaterialEntry::get_material_id() const { return material_id; }

} // namespace godot
