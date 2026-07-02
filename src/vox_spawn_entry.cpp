#include "vox_spawn_entry.h"

#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/core/math.hpp>

namespace godot {

VoxSpawnEntry::VoxSpawnEntry() {}
VoxSpawnEntry::~VoxSpawnEntry() {}

void VoxSpawnEntry::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_palette_index", "index"), &VoxSpawnEntry::set_palette_index);
    ClassDB::bind_method(D_METHOD("get_palette_index"), &VoxSpawnEntry::get_palette_index);
    ClassDB::bind_method(D_METHOD("set_spawn_type", "type"), &VoxSpawnEntry::set_spawn_type);
    ClassDB::bind_method(D_METHOD("get_spawn_type"), &VoxSpawnEntry::get_spawn_type);

    ADD_PROPERTY(PropertyInfo(Variant::INT, "palette_index", PROPERTY_HINT_RANGE, "0,255,1"), "set_palette_index", "get_palette_index");
    ADD_PROPERTY(PropertyInfo(Variant::STRING, "spawn_type"), "set_spawn_type", "get_spawn_type");
}

void VoxSpawnEntry::set_palette_index(int p_index) {
    palette_index = Math::clamp(p_index, 0, 255);
    emit_changed();
}

int VoxSpawnEntry::get_palette_index() const {
    return palette_index;
}

void VoxSpawnEntry::set_spawn_type(const String &p_type) {
    spawn_type = p_type;
    emit_changed();
}

String VoxSpawnEntry::get_spawn_type() const {
    return spawn_type;
}

} // namespace godot
