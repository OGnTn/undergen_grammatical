#ifndef VOX_SPAWN_ENTRY_H
#define VOX_SPAWN_ENTRY_H

#include <godot_cpp/classes/resource.hpp>

namespace godot {

class VoxSpawnEntry : public Resource {
    GDCLASS(VoxSpawnEntry, Resource);

private:
    int palette_index = 0;
    String spawn_type;

protected:
    static void _bind_methods();

public:
    VoxSpawnEntry();
    virtual ~VoxSpawnEntry();

    void set_palette_index(int p_index);
    int get_palette_index() const;

    void set_spawn_type(const String &p_type);
    String get_spawn_type() const;
};

} // namespace godot

#endif // VOX_SPAWN_ENTRY_H
