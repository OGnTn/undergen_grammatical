#ifndef UNDERGEN_MESH_STAMP_NODE_H
#define UNDERGEN_MESH_STAMP_NODE_H

#include "undergen_node.h"
#include "density_grid.h"
#include "level_gen_data.h"

#include <godot_cpp/classes/mesh.hpp>
#include <godot_cpp/classes/array_mesh.hpp>
#include <godot_cpp/classes/resource_loader.hpp>
#include <godot_cpp/variant/dictionary.hpp>
#include <godot_cpp/variant/transform3d.hpp>
#include <godot_cpp/variant/vector3.hpp>
#include <godot_cpp/variant/quaternion.hpp>
#include <godot_cpp/variant/basis.hpp>

#include <vector>
#include <map>

namespace godot {

struct MeshStampTriangle {
    Vector3 a, b, c;
    Vector3 na, nb, nc;
    Vector3 face_normal;
    bool is_sharp = false;
    Vector3 min_bounds;
    Vector3 max_bounds;
};

class UnderGenMeshStampNode : public UnderGenNode {
    GDCLASS(UnderGenMeshStampNode, UnderGenNode);

public:
    enum BlendMode {
        MODE_OVERWRITE = 0,
        MODE_MAX_UNION = 1,
        MODE_SUBTRACT_CARVE = 2
    };

private:
    Ref<Mesh> mesh;
    String mesh_path;
    int material_id = 1;
    int blend_mode = MODE_OVERWRITE;
    float sharp_edge_angle_threshold = 45.0f; // degrees
    bool force_sharp_normals = false;

    Vector3 position_offset = Vector3(0, 0, 0);
    Vector3 rotation_degrees = Vector3(0, 0, 0);
    Vector3 scale = Vector3(1, 1, 1);

    bool exclude_from_smoothing = false;
    bool exclude_from_warping = false;
    bool clear_room_air = true;

    // Internal cache for loaded mesh path
    Ref<Mesh> _cached_loaded_mesh;
    String _last_loaded_path;

    Ref<Mesh> _get_effective_mesh();

protected:
    static void _bind_methods();

public:
    UnderGenMeshStampNode();
    virtual ~UnderGenMeshStampNode();

    static void stamp_mesh_resource(
        DensityGrid *grid,
        const Ref<Mesh> &p_mesh,
        const Transform3D &p_transform,
        int p_material_id,
        int p_blend_mode,
        float p_sharp_threshold_deg,
        bool p_force_sharp_normals,
        int p_zone_id = 0,
        bool p_clear_air = false,
        Vector3i p_room_origin = Vector3i(),
        Vector3i p_room_size = Vector3i()
    );

    void set_mesh(const Ref<Mesh> &p_mesh);
    Ref<Mesh> get_mesh() const;

    void set_mesh_path(const String &p_path);
    String get_mesh_path() const;

    void set_material_id(int p_id);
    int get_material_id() const;

    void set_blend_mode(int p_mode);
    int get_blend_mode() const;

    void set_sharp_edge_angle_threshold(float p_angle_deg);
    float get_sharp_edge_angle_threshold() const;

    void set_force_sharp_normals(bool p_force);
    bool get_force_sharp_normals() const;

    void set_position_offset(const Vector3 &p_offset);
    Vector3 get_position_offset() const;

    void set_rotation_degrees(const Vector3 &p_rot);
    Vector3 get_rotation_degrees() const;

    void set_scale(const Vector3 &p_scale);
    Vector3 get_scale() const;

    void set_exclude_from_smoothing(bool p_exclude);
    bool get_exclude_from_smoothing() const;

    void set_exclude_from_warping(bool p_exclude);
    bool get_exclude_from_warping() const;

    void set_clear_room_air(bool p_clear);
    bool get_clear_room_air() const;

    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
};

} // namespace godot

VARIANT_ENUM_CAST(UnderGenMeshStampNode::BlendMode);

#endif // UNDERGEN_MESH_STAMP_NODE_H
