#ifndef ZONE_MATERIAL_ENTRY_H
#define ZONE_MATERIAL_ENTRY_H

#include <godot_cpp/classes/resource.hpp>

namespace godot {

// Simple Resource holding a single zone-name → material-ID mapping.
// Used in TypedArrays so Godot's inspector renders editable list entries
// instead of the broken untyped Dictionary editor.
class ZoneMaterialEntry : public Resource {
    GDCLASS(ZoneMaterialEntry, Resource);

private:
    String zone_name;
    int material_id = 0;

protected:
    static void _bind_methods();

public:
    ZoneMaterialEntry();
    virtual ~ZoneMaterialEntry();

    void set_zone_name(const String &p_name);
    String get_zone_name() const;

    void set_material_id(int p_id);
    int get_material_id() const;
};

} // namespace godot

#endif // ZONE_MATERIAL_ENTRY_H
