#ifndef VOX_MATERIAL_ENTRY_H
#define VOX_MATERIAL_ENTRY_H

#include <godot_cpp/classes/resource.hpp>

namespace godot {

class VoxMaterialEntry : public Resource {
    GDCLASS(VoxMaterialEntry, Resource);

private:
    int palette_index = 0;
    int material_id = 0;
    float thickness = 1.0f;
    bool stepped = true;

protected:
    static void _bind_methods();

public:
    VoxMaterialEntry();
    virtual ~VoxMaterialEntry();

    void set_palette_index(int p_index);
    int get_palette_index() const;

    void set_material_id(int p_id);
    int get_material_id() const;

    void set_thickness(float p_thickness);
    float get_thickness() const;

    void set_stepped(bool p_stepped);
    bool is_stepped() const;
};

} // namespace godot

#endif // VOX_MATERIAL_ENTRY_H
