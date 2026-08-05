#include "undergen_mesh_stamp_node.h"
#include <godot_cpp/classes/file_access.hpp>
#include <godot_cpp/variant/utility_functions.hpp>
#include <godot_cpp/core/math.hpp>
#include <algorithm>
#include <cmath>
#include <cfloat>

namespace godot {

// --- Helper Functions for Geometry & Intersection ---

static Vector3 closest_point_on_triangle(const Vector3 &p, const Vector3 &a, const Vector3 &b, const Vector3 &c) {
    Vector3 ab = b - a;
    Vector3 ac = c - a;
    Vector3 ap = p - a;

    float d1 = ab.dot(ap);
    float d2 = ac.dot(ap);
    if (d1 <= 0.0f && d2 <= 0.0f) return a;

    Vector3 bp = p - b;
    float d3 = ab.dot(bp);
    float d4 = ac.dot(bp);
    if (d3 >= 0.0f && d4 <= d3) return b;

    float vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {
        float v = d1 / (d1 - d3);
        return a + v * ab;
    }

    Vector3 cp = p - c;
    float d5 = ab.dot(cp);
    float d6 = ac.dot(cp);
    if (d6 >= 0.0f && d5 <= d6) return c;

    float vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
        float w = d2 / (d2 - d6);
        return a + w * ac;
    }

    float va = d3 * d6 - d5 * d4;
    if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
        float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + w * (c - b);
    }

    float denom = 1.0f / (va + vb + vc);
    float v = vb * denom;
    float w = vc * denom;
    return a + ab * v + ac * w;
}

static bool ray_triangle_intersect(
    const Vector3 &ray_origin, const Vector3 &ray_dir,
    const Vector3 &v0, const Vector3 &v1, const Vector3 &v2,
    float &out_t, float &out_u, float &out_v)
{
    const float EPSILON = 1e-7f;
    Vector3 edge1 = v1 - v0;
    Vector3 edge2 = v2 - v0;
    Vector3 h = ray_dir.cross(edge2);
    float a = edge1.dot(h);

    if (a > -EPSILON && a < EPSILON) return false; // Ray is parallel to triangle

    float f = 1.0f / a;
    Vector3 s = ray_origin - v0;
    out_u = f * s.dot(h);
    if (out_u < 0.0f || out_u > 1.0f) return false;

    Vector3 q = s.cross(edge1);
    out_v = f * ray_dir.dot(q);
    if (out_v < 0.0f || out_u + out_v > 1.0f) return false;

    out_t = f * edge2.dot(q);
    return out_t > EPSILON;
}

static bool segment_triangle_intersect(
    const Vector3 &p0, const Vector3 &p1,
    const Vector3 &v0, const Vector3 &v1, const Vector3 &v2,
    float &out_t, float &out_u, float &out_v)
{
    Vector3 dir = p1 - p0;
    float len = dir.length();
    if (len < 1e-7f) return false;

    Vector3 ray_dir = dir / len;
    float t_hit = 0.0f;
    if (ray_triangle_intersect(p0, ray_dir, v0, v1, v2, t_hit, out_u, out_v)) {
        if (t_hit <= len) {
            out_t = t_hit / len; // Normalize t to [0, 1]
            return true;
        }
    }
    return false;
}

UnderGenMeshStampNode::UnderGenMeshStampNode() {}

UnderGenMeshStampNode::~UnderGenMeshStampNode() {}

void UnderGenMeshStampNode::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_mesh", "mesh"), &UnderGenMeshStampNode::set_mesh);
    ClassDB::bind_method(D_METHOD("get_mesh"), &UnderGenMeshStampNode::get_mesh);

    ClassDB::bind_method(D_METHOD("set_mesh_path", "path"), &UnderGenMeshStampNode::set_mesh_path);
    ClassDB::bind_method(D_METHOD("get_mesh_path"), &UnderGenMeshStampNode::get_mesh_path);

    ClassDB::bind_method(D_METHOD("set_material_id", "id"), &UnderGenMeshStampNode::set_material_id);
    ClassDB::bind_method(D_METHOD("get_material_id"), &UnderGenMeshStampNode::get_material_id);

    ClassDB::bind_method(D_METHOD("set_blend_mode", "mode"), &UnderGenMeshStampNode::set_blend_mode);
    ClassDB::bind_method(D_METHOD("get_blend_mode"), &UnderGenMeshStampNode::get_blend_mode);

    ClassDB::bind_method(D_METHOD("set_sharp_edge_angle_threshold", "angle"), &UnderGenMeshStampNode::set_sharp_edge_angle_threshold);
    ClassDB::bind_method(D_METHOD("get_sharp_edge_angle_threshold"), &UnderGenMeshStampNode::get_sharp_edge_angle_threshold);

    ClassDB::bind_method(D_METHOD("set_force_sharp_normals", "force"), &UnderGenMeshStampNode::set_force_sharp_normals);
    ClassDB::bind_method(D_METHOD("get_force_sharp_normals"), &UnderGenMeshStampNode::get_force_sharp_normals);

    ClassDB::bind_method(D_METHOD("set_position_offset", "offset"), &UnderGenMeshStampNode::set_position_offset);
    ClassDB::bind_method(D_METHOD("get_position_offset"), &UnderGenMeshStampNode::get_position_offset);

    ClassDB::bind_method(D_METHOD("set_rotation_degrees", "rot"), &UnderGenMeshStampNode::set_rotation_degrees);
    ClassDB::bind_method(D_METHOD("get_rotation_degrees"), &UnderGenMeshStampNode::get_rotation_degrees);

    ClassDB::bind_method(D_METHOD("set_scale", "scale"), &UnderGenMeshStampNode::set_scale);
    ClassDB::bind_method(D_METHOD("get_scale"), &UnderGenMeshStampNode::get_scale);

    ClassDB::bind_method(D_METHOD("set_exclude_from_smoothing", "exclude"), &UnderGenMeshStampNode::set_exclude_from_smoothing);
    ClassDB::bind_method(D_METHOD("get_exclude_from_smoothing"), &UnderGenMeshStampNode::get_exclude_from_smoothing);

    ClassDB::bind_method(D_METHOD("set_exclude_from_warping", "exclude"), &UnderGenMeshStampNode::set_exclude_from_warping);
    ClassDB::bind_method(D_METHOD("get_exclude_from_warping"), &UnderGenMeshStampNode::get_exclude_from_warping);

    ClassDB::bind_method(D_METHOD("set_clear_room_air", "clear"), &UnderGenMeshStampNode::set_clear_room_air);
    ClassDB::bind_method(D_METHOD("get_clear_room_air"), &UnderGenMeshStampNode::get_clear_room_air);

    ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "mesh", PROPERTY_HINT_RESOURCE_TYPE, "Mesh"), "set_mesh", "get_mesh");
    ADD_PROPERTY(PropertyInfo(Variant::STRING, "mesh_path", PROPERTY_HINT_FILE, "*.res,*.tres,*.obj,*.gltf,*.glb"), "set_mesh_path", "get_mesh_path");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "material_id"), "set_material_id", "get_material_id");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "blend_mode", PROPERTY_HINT_ENUM, "Overwrite,MaxUnion,SubtractCarve"), "set_blend_mode", "get_blend_mode");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "sharp_edge_angle_threshold", PROPERTY_HINT_RANGE, "0.0, 180.0, 0.5"), "set_sharp_edge_angle_threshold", "get_sharp_edge_angle_threshold");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "force_sharp_normals"), "set_force_sharp_normals", "get_force_sharp_normals");

    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3, "position_offset"), "set_position_offset", "get_position_offset");
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3, "rotation_degrees"), "set_rotation_degrees", "get_rotation_degrees");
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3, "scale"), "set_scale", "get_scale");

    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "exclude_from_smoothing"), "set_exclude_from_smoothing", "get_exclude_from_smoothing");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "exclude_from_warping"), "set_exclude_from_warping", "get_exclude_from_warping");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "clear_room_air"), "set_clear_room_air", "get_clear_room_air");

    BIND_ENUM_CONSTANT(MODE_OVERWRITE);
    BIND_ENUM_CONSTANT(MODE_MAX_UNION);
    BIND_ENUM_CONSTANT(MODE_SUBTRACT_CARVE);
}

void UnderGenMeshStampNode::set_mesh(const Ref<Mesh> &p_mesh) { mesh = p_mesh; }
Ref<Mesh> UnderGenMeshStampNode::get_mesh() const { return mesh; }

void UnderGenMeshStampNode::set_mesh_path(const String &p_path) { mesh_path = p_path; }
String UnderGenMeshStampNode::get_mesh_path() const { return mesh_path; }

void UnderGenMeshStampNode::set_material_id(int p_id) { material_id = p_id; }
int UnderGenMeshStampNode::get_material_id() const { return material_id; }

void UnderGenMeshStampNode::set_blend_mode(int p_mode) { blend_mode = p_mode; }
int UnderGenMeshStampNode::get_blend_mode() const { return blend_mode; }

void UnderGenMeshStampNode::set_sharp_edge_angle_threshold(float p_angle_deg) { sharp_edge_angle_threshold = p_angle_deg; }
float UnderGenMeshStampNode::get_sharp_edge_angle_threshold() const { return sharp_edge_angle_threshold; }

void UnderGenMeshStampNode::set_force_sharp_normals(bool p_force) { force_sharp_normals = p_force; }
bool UnderGenMeshStampNode::get_force_sharp_normals() const { return force_sharp_normals; }

void UnderGenMeshStampNode::set_position_offset(const Vector3 &p_offset) { position_offset = p_offset; }
Vector3 UnderGenMeshStampNode::get_position_offset() const { return position_offset; }

void UnderGenMeshStampNode::set_rotation_degrees(const Vector3 &p_rot) { rotation_degrees = p_rot; }
Vector3 UnderGenMeshStampNode::get_rotation_degrees() const { return rotation_degrees; }

void UnderGenMeshStampNode::set_scale(const Vector3 &p_scale) { scale = p_scale; }
Vector3 UnderGenMeshStampNode::get_scale() const { return scale; }

void UnderGenMeshStampNode::set_exclude_from_smoothing(bool p_exclude) { exclude_from_smoothing = p_exclude; }
bool UnderGenMeshStampNode::get_exclude_from_smoothing() const { return exclude_from_smoothing; }

void UnderGenMeshStampNode::set_exclude_from_warping(bool p_exclude) { exclude_from_warping = p_exclude; }
bool UnderGenMeshStampNode::get_exclude_from_warping() const { return exclude_from_warping; }

void UnderGenMeshStampNode::set_clear_room_air(bool p_clear) { clear_room_air = p_clear; }
bool UnderGenMeshStampNode::get_clear_room_air() const { return clear_room_air; }

Ref<Mesh> UnderGenMeshStampNode::_get_effective_mesh() {
    if (mesh.is_valid()) return mesh;
    if (!mesh_path.is_empty()) {
        if (_last_loaded_path == mesh_path && _cached_loaded_mesh.is_valid()) {
            return _cached_loaded_mesh;
        }
        ResourceLoader *loader = ResourceLoader::get_singleton();
        if (loader && loader->exists(mesh_path)) {
            _cached_loaded_mesh = loader->load(mesh_path);
            _last_loaded_path = mesh_path;
            return _cached_loaded_mesh;
        }
    }
    return Ref<Mesh>();
}

void UnderGenMeshStampNode::stamp_mesh_resource(
    DensityGrid *grid,
    const Ref<Mesh> &p_mesh,
    const Transform3D &p_transform,
    int p_material_id,
    int p_blend_mode,
    float p_sharp_threshold_deg,
    bool p_force_sharp_normals,
    int p_zone_id,
    bool p_clear_air,
    Vector3i p_room_origin,
    Vector3i p_room_size)
{
    if (!grid || p_mesh.is_null()) return;

    // 1. Extract Triangles from Godot Mesh
    std::vector<MeshStampTriangle> triangles;
    int surf_count = p_mesh->get_surface_count();

    for (int s = 0; s < surf_count; ++s) {
        Array arrays = p_mesh->surface_get_arrays(s);
        if (arrays.size() <= Mesh::ARRAY_VERTEX) continue;

        PackedVector3Array verts = arrays[Mesh::ARRAY_VERTEX];
        if (verts.size() < 3) continue;

        PackedVector3Array norms;
        if (arrays.size() > Mesh::ARRAY_NORMAL && arrays[Mesh::ARRAY_NORMAL].get_type() == Variant::PACKED_VECTOR3_ARRAY) {
            norms = arrays[Mesh::ARRAY_NORMAL];
        }

        PackedInt32Array indices;
        if (arrays.size() > Mesh::ARRAY_INDEX && arrays[Mesh::ARRAY_INDEX].get_type() == Variant::PACKED_INT32_ARRAY) {
            indices = arrays[Mesh::ARRAY_INDEX];
        }

        bool has_indices = indices.size() > 0;
        int num_triangles = has_indices ? (indices.size() / 3) : (verts.size() / 3);

        Basis normal_basis = p_transform.basis.inverse().transposed();

        for (int i = 0; i < num_triangles; ++i) {
            int idx0 = has_indices ? indices[i * 3 + 0] : (i * 3 + 0);
            int idx1 = has_indices ? indices[i * 3 + 1] : (i * 3 + 1);
            int idx2 = has_indices ? indices[i * 3 + 2] : (i * 3 + 2);

            MeshStampTriangle tri;
            tri.a = p_transform.xform(verts[idx0]);
            tri.b = p_transform.xform(verts[idx1]);
            tri.c = p_transform.xform(verts[idx2]);

            Vector3 edge1 = tri.b - tri.a;
            Vector3 edge2 = tri.c - tri.a;
            Vector3 fn = edge1.cross(edge2);
            if (fn.length_squared() < 1e-12f) continue; // Degenerate triangle
            tri.face_normal = fn.normalized();

            if (norms.size() > idx2) {
                tri.na = normal_basis.xform(norms[idx0]).normalized();
                tri.nb = normal_basis.xform(norms[idx1]).normalized();
                tri.nc = normal_basis.xform(norms[idx2]).normalized();
            } else {
                tri.na = tri.face_normal;
                tri.nb = tri.face_normal;
                tri.nc = tri.face_normal;
            }

            // Sharpness check
            if (p_force_sharp_normals) {
                tri.is_sharp = true;
            } else {
                float cos_thresh = std::cos(Math::deg_to_rad(p_sharp_threshold_deg));
                float dot_a = std::abs(tri.na.dot(tri.face_normal));
                float dot_b = std::abs(tri.nb.dot(tri.face_normal));
                float dot_c = std::abs(tri.nc.dot(tri.face_normal));
                if (dot_a < cos_thresh || dot_b < cos_thresh || dot_c < cos_thresh) {
                    tri.is_sharp = true;
                }
            }

            tri.min_bounds = Vector3(
                std::min({tri.a.x, tri.b.x, tri.c.x}),
                std::min({tri.a.y, tri.b.y, tri.c.y}),
                std::min({tri.a.z, tri.b.z, tri.c.z})
            );
            tri.max_bounds = Vector3(
                std::max({tri.a.x, tri.b.x, tri.c.x}),
                std::max({tri.a.y, tri.b.y, tri.c.y}),
                std::max({tri.a.z, tri.b.z, tri.c.z})
            );

            triangles.push_back(tri);
        }
    }

    if (triangles.empty()) return;

    // 2. Clear room air if requested
    if (p_clear_air && p_room_size.x > 0 && p_room_size.y > 0 && p_room_size.z > 0) {
        for (int rx = 0; rx < p_room_size.x; ++rx) {
            for (int ry = 0; ry < p_room_size.y; ++ry) {
                for (int rz = 0; rz < p_room_size.z; ++rz) {
                    Vector3i cpos = p_room_origin + Vector3i(rx, ry, rz);
                    if (grid->is_valid_position(cpos)) {
                        grid->set_cell(cpos, 0.0f);
                        if (p_zone_id > 0) grid->set_zone_at(cpos, p_zone_id);
                    }
                }
            }
        }
    }

    // 3. Compute AABB of all transformed mesh triangles
    Vector3 mesh_aabb_min = triangles[0].min_bounds;
    Vector3 mesh_aabb_max = triangles[0].max_bounds;
    for (size_t i = 1; i < triangles.size(); ++i) {
        mesh_aabb_min.x = std::min(mesh_aabb_min.x, triangles[i].min_bounds.x);
        mesh_aabb_min.y = std::min(mesh_aabb_min.y, triangles[i].min_bounds.y);
        mesh_aabb_min.z = std::min(mesh_aabb_min.z, triangles[i].min_bounds.z);

        mesh_aabb_max.x = std::max(mesh_aabb_max.x, triangles[i].max_bounds.x);
        mesh_aabb_max.y = std::max(mesh_aabb_max.y, triangles[i].max_bounds.y);
        mesh_aabb_max.z = std::max(mesh_aabb_max.z, triangles[i].max_bounds.z);
    }

    // Expand bounding box slightly for voxel neighborhood
    int min_gx = (int)std::floor(mesh_aabb_min.x) - 2;
    int max_gx = (int)std::ceil(mesh_aabb_max.x) + 2;
    int min_gy = (int)std::floor(mesh_aabb_min.y) - 2;
    int max_gy = (int)std::ceil(mesh_aabb_max.y) + 2;
    int min_gz = (int)std::floor(mesh_aabb_min.z) - 2;
    int max_gz = (int)std::ceil(mesh_aabb_max.z) + 2;

    Vector3i gdim = grid->get_grid_dimensions();
    min_gx = std::max(0, min_gx); max_gx = std::min(gdim.x - 1, max_gx);
    min_gy = std::max(0, min_gy); max_gy = std::min(gdim.y - 1, max_gy);
    min_gz = std::max(0, min_gz); max_gz = std::min(gdim.z - 1, max_gz);

    if (min_gx > max_gx || min_gy > max_gy || min_gz > max_gz) return;

    // 4. Voxel Density Computation (Ray Parity + Point-Triangle Distance)
    std::unordered_map<Vector3i, bool, Vector3iHash> solid_mask;
    std::unordered_map<Vector3i, float, Vector3iHash> density_map;

    Vector3 ray_dir(1.00001f, 0.00003f, 0.00007f); // Slightly off-axis ray to avoid parallel edge degeneracy
    ray_dir.normalize();

    for (int gx = min_gx; gx <= max_gx; ++gx) {
        for (int gy = min_gy; gy <= max_gy; ++gy) {
            for (int gz = min_gz; gz <= max_gz; ++gz) {
                Vector3i gpos(gx, gy, gz);
                Vector3 p(gpos);

                // Raycast to determine inside/outside
                int hit_count = 0;
                float min_dist_sq = FLT_MAX;

                for (const auto &tri : triangles) {
                    float t_hit = 0, u = 0, v = 0;
                    if (ray_triangle_intersect(p, ray_dir, tri.a, tri.b, tri.c, t_hit, u, v)) {
                        hit_count++;
                    }

                    Vector3 cp = closest_point_on_triangle(p, tri.a, tri.b, tri.c);
                    float dsq = p.distance_squared_to(cp);
                    if (dsq < min_dist_sq) {
                        min_dist_sq = dsq;
                    }
                }

                bool is_inside = (hit_count % 2) == 1;
                float min_dist = std::sqrt(min_dist_sq);
                float sdf = is_inside ? -min_dist : min_dist;

                // Density mapping: > 0.0 is solid (inside), <= 0.0 is air (outside)
                // At sdf = 0.0 (on surface), density = 0.0 (exact surface threshold)
                float density = Math::clamp(-sdf, -1.0f, 1.0f);

                solid_mask[gpos] = is_inside || (density > 0.0f);
                density_map[gpos] = density;

                if (p_blend_mode == MODE_OVERWRITE) {
                    grid->set_cell(gpos, density);
                } else if (p_blend_mode == MODE_MAX_UNION) {
                    float cur = grid->get_cell(gpos);
                    grid->set_cell(gpos, std::max(cur, density));
                } else if (p_blend_mode == MODE_SUBTRACT_CARVE) {
                    float cur = grid->get_cell(gpos);
                    grid->set_cell(gpos, std::min(cur, 1.0f - density));
                }

                if (is_inside || density > 0.0f) {
                    grid->set_material_id(gpos, p_material_id);
                }
                if (p_zone_id > 0) {
                    grid->set_zone_at(gpos, p_zone_id);
                }
            }
        }
    }

    // 5. Hermite Edge Intersection & Normal Computation (For Dual Contouring Sharp & Smooth Edges!)
    const Vector3i axes[3] = { Vector3i(1, 0, 0), Vector3i(0, 1, 0), Vector3i(0, 0, 1) };

    for (int gx = min_gx; gx <= max_gx; ++gx) {
        for (int gy = min_gy; gy <= max_gy; ++gy) {
            for (int gz = min_gz; gz <= max_gz; ++gz) {
                Vector3i p1_pos(gx, gy, gz);
                bool s1 = solid_mask[p1_pos];

                for (int axis = 0; axis < 3; ++axis) {
                    Vector3i p2_pos = p1_pos + axes[axis];
                    if (p2_pos.x > max_gx || p2_pos.y > max_gy || p2_pos.z > max_gz) continue;

                    bool s2 = solid_mask[p2_pos];
                    if (s1 != s2) { // Isosurface crossing edge
                        Vector3 p1(p1_pos);
                        Vector3 p2(p2_pos);

                        float best_t = 2.0f;
                        Vector3 best_normal(0, 1, 0);

                        for (const auto &tri : triangles) {
                            float t = 0, u = 0, v = 0;
                            if (segment_triangle_intersect(p1, p2, tri.a, tri.b, tri.c, t, u, v)) {
                                if (t >= 0.0f && t <= 1.0f && t < best_t) {
                                    best_t = t;
                                    float w = 1.0f - u - v;

                                    Vector3 n_smooth = (w * tri.na + u * tri.nb + v * tri.nc);
                                    if (n_smooth.length_squared() > 1e-8f) {
                                        n_smooth.normalize();
                                    } else {
                                        n_smooth = tri.face_normal;
                                    }

                                    if (tri.is_sharp || p_force_sharp_normals) {
                                        best_normal = tri.face_normal;
                                    } else {
                                        best_normal = n_smooth;
                                    }
                                }
                            }
                        }

                        if (best_t <= 1.0f) {
                            grid->set_hermite_edge(p1_pos, axis, best_t, best_normal);
                        }
                    }
                }
            }
        }
    }
}

void UnderGenMeshStampNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    Dictionary context = inputs.get(0, Dictionary());
    if (context.is_empty()) return;

    Ref<DensityGrid> grid = context.get("grid", Ref<DensityGrid>());
    if (grid.is_null()) { outputs[0] = context; return; }

    Ref<Mesh> effective_mesh = _get_effective_mesh();
    if (effective_mesh.is_null()) {
        UtilityFunctions::printerr("UnderGenMeshStampNode: No valid mesh or mesh_path provided.");
        outputs[0] = context;
        return;
    }

    Transform3D xform;
    xform.origin = position_offset;
    Vector3 rot_rad = Vector3(
        Math::deg_to_rad(rotation_degrees.x),
        Math::deg_to_rad(rotation_degrees.y),
        Math::deg_to_rad(rotation_degrees.z)
    );
    xform.basis = Basis::from_euler(rot_rad).scaled(scale);

    Array rooms_arr = context.get("rooms", Array());
    if (rooms_arr.size() > 0) {
        for (int i = 0; i < rooms_arr.size(); ++i) {
            Dictionary r_dict = rooms_arr[i];
            ResolvedRoom room;
            room.id = r_dict.get("id", "");
            room.type = r_dict.get("type", "");
            room.position = r_dict.get("position", Vector3i());
            room.size = r_dict.get("size", Vector3i());

            Transform3D room_xform = xform;
            room_xform.origin += Vector3(room.position);

            int zone_id = grid->register_zone_name(room.type);
            stamp_mesh_resource(
                grid.ptr(), effective_mesh, room_xform,
                material_id, blend_mode, sharp_edge_angle_threshold, force_sharp_normals,
                zone_id, clear_room_air, room.position, room.size
            );

            if (exclude_from_smoothing) r_dict["exclude_from_smoothing"] = true;
            if (exclude_from_warping) r_dict["exclude_from_warping"] = true;
        }
    } else {
        stamp_mesh_resource(
            grid.ptr(), effective_mesh, xform,
            material_id, blend_mode, sharp_edge_angle_threshold, force_sharp_normals
        );
    }

    outputs[0] = context;
}

} // namespace godot
