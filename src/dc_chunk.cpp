// src/dc_chunk.cpp
#include "dc_chunk.h"
#include "mc_tables.h"

#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/variant/utility_functions.hpp>
#include <godot_cpp/core/math.hpp>
#include <godot_cpp/classes/node.hpp>
#include <godot_cpp/classes/static_body3d.hpp>
#include <godot_cpp/classes/array_mesh.hpp>
#include <godot_cpp/classes/occluder_instance3d.hpp>
#include <godot_cpp/classes/array_occluder3d.hpp>
#include <godot_cpp/classes/time.hpp>
#include <godot_cpp/classes/area3d.hpp>

#include <unordered_map>
#include <vector>
#include <map>

namespace godot {

// Simple hash function helper for Vector3i keys
struct Vector3iHash {
    std::size_t operator()(const Vector3i& v) const {
        return (std::hash<int>()(v.x) ^ (std::hash<int>()(v.y) << 1)) ^ (std::hash<int>()(v.z) << 2);
    }
};

struct CellVertexInfo {
    Vector3 position;
    Vector3 normal;
    int material_id = 0;
    int master_index = -1;
};

struct SurfaceBuilder {
    PackedVector3Array vertices;
    PackedVector3Array normals;
    PackedInt32Array indices;
    std::map<int, int> index_cache;
};

DCChunk::DCChunk() {
    Ref<ArrayMesh> new_mesh;
    new_mesh.instantiate();
    set_mesh(new_mesh);
}

DCChunk::~DCChunk() {
    _clear_collision();
    _clear_occluder();
    _clear_liquid();
}

void DCChunk::_notification(int p_what) {}

void DCChunk::_bind_methods() {
    ClassDB::bind_method(D_METHOD("generate_mesh_from_density_grid"), &DCChunk::generate_mesh_from_density_grid);

    ClassDB::bind_method(D_METHOD("set_chunk_size", "size"), &DCChunk::set_chunk_size);
    ClassDB::bind_method(D_METHOD("get_chunk_size"), &DCChunk::get_chunk_size);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "chunk_size"), "set_chunk_size", "get_chunk_size");

    ClassDB::bind_method(D_METHOD("set_voxel_size", "size"), &DCChunk::set_voxel_size);
    ClassDB::bind_method(D_METHOD("get_voxel_size"), &DCChunk::get_voxel_size);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "voxel_size"), "set_voxel_size", "get_voxel_size");

    ClassDB::bind_method(D_METHOD("set_chunk_grid_offset", "offset"), &DCChunk::set_chunk_grid_offset);
    ClassDB::bind_method(D_METHOD("get_chunk_grid_offset"), &DCChunk::get_chunk_grid_offset);
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3I, "chunk_grid_offset"), "set_chunk_grid_offset", "get_chunk_grid_offset");

    ClassDB::bind_method(D_METHOD("set_materials", "materials"), &DCChunk::set_materials);
    ClassDB::bind_method(D_METHOD("get_materials"), &DCChunk::get_materials);
    ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "materials", PROPERTY_HINT_TYPE_STRING, "4/17:Material"), "set_materials", "get_materials");

    ClassDB::bind_method(D_METHOD("set_density_grid", "grid"), &DCChunk::set_density_grid);
    ClassDB::bind_method(D_METHOD("get_density_grid"), &DCChunk::get_density_grid);
    ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "density_grid", PROPERTY_HINT_RESOURCE_TYPE, "DensityGrid"), "set_density_grid", "get_density_grid");

    ClassDB::bind_method(D_METHOD("set_generate_collision", "enable"), &DCChunk::set_generate_collision);
    ClassDB::bind_method(D_METHOD("get_generate_collision"), &DCChunk::get_generate_collision);
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "generate_collision"), "set_generate_collision", "get_generate_collision");

    ClassDB::bind_method(D_METHOD("set_generate_occluder", "enable"), &DCChunk::set_generate_occluder);
    ClassDB::bind_method(D_METHOD("get_generate_occluder"), &DCChunk::get_generate_occluder);
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "generate_occluder"), "set_generate_occluder", "get_generate_occluder");

    ClassDB::bind_method(D_METHOD("set_compute_shader", "shader"), &DCChunk::set_compute_shader);
    ClassDB::bind_method(D_METHOD("get_compute_shader"), &DCChunk::get_compute_shader);
    ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "compute_shader", PROPERTY_HINT_RESOURCE_TYPE, "RDShaderFile"), "set_compute_shader", "get_compute_shader");

    ClassDB::bind_method(D_METHOD("set_liquid_material", "material"), &DCChunk::set_liquid_material);
    ClassDB::bind_method(D_METHOD("get_liquid_material"), &DCChunk::get_liquid_material);
    ClassDB::bind_method(D_METHOD("set_liquid_material_id", "id"), &DCChunk::set_liquid_material_id);
    ClassDB::bind_method(D_METHOD("get_liquid_material_id"), &DCChunk::get_liquid_material_id);
    ClassDB::bind_method(D_METHOD("set_generate_liquid_trigger", "enable"), &DCChunk::set_generate_liquid_trigger);
    ClassDB::bind_method(D_METHOD("get_generate_liquid_trigger"), &DCChunk::get_generate_liquid_trigger);
    ClassDB::bind_method(D_METHOD("set_flow_spread_limit", "limit"), &DCChunk::set_flow_spread_limit);
    ClassDB::bind_method(D_METHOD("get_flow_spread_limit"), &DCChunk::get_flow_spread_limit);

    ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "liquid_material", PROPERTY_HINT_RESOURCE_TYPE, "Material"), "set_liquid_material", "get_liquid_material");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "liquid_material_id", PROPERTY_HINT_RANGE, "0,255,1"), "set_liquid_material_id", "get_liquid_material_id");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "generate_liquid_trigger"), "set_generate_liquid_trigger", "get_generate_liquid_trigger");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "flow_spread_limit"), "set_flow_spread_limit", "get_flow_spread_limit");

    ClassDB::bind_method(D_METHOD("set_smooth_normals", "enable"), &DCChunk::set_smooth_normals);
    ClassDB::bind_method(D_METHOD("get_smooth_normals"), &DCChunk::get_smooth_normals);
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "smooth_normals"), "set_smooth_normals", "get_smooth_normals");

    ClassDB::bind_method(D_METHOD("set_flip_normals", "enable"), &DCChunk::set_flip_normals);
    ClassDB::bind_method(D_METHOD("get_flip_normals"), &DCChunk::get_flip_normals);
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "flip_normals"), "set_flip_normals", "get_flip_normals");

    // Dual Contouring Specific
    ClassDB::bind_method(D_METHOD("set_use_qef", "use"), &DCChunk::set_use_qef);
    ClassDB::bind_method(D_METHOD("get_use_qef"), &DCChunk::get_use_qef);
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "use_qef"), "set_use_qef", "get_use_qef");

    ClassDB::bind_method(D_METHOD("set_qef_regularization", "regularization"), &DCChunk::set_qef_regularization);
    ClassDB::bind_method(D_METHOD("get_qef_regularization"), &DCChunk::get_qef_regularization);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "qef_regularization"), "set_qef_regularization", "get_qef_regularization");
}

Vector3 solve_qef_cramer(const std::vector<Vector3> &pts, const std::vector<Vector3> &norms, const Vector3 &cell_min, const Vector3 &cell_max, float regularization) {
    int k = pts.size();
    if (k == 0) return (cell_min + cell_max) * 0.5f;

    Vector3 fallback(0, 0, 0);
    for (const Vector3 &p : pts) fallback += p;
    fallback /= (float)k;

    double mxx = 0, mxy = 0, mxz = 0;
    double myy = 0, myz = 0, mzz = 0;
    double yx = 0, yy = 0, yz = 0;

    for (int i = 0; i < k; ++i) {
        double nx = norms[i].x;
        double ny = norms[i].y;
        double nz = norms[i].z;
        double dot = nx * pts[i].x + ny * pts[i].y + nz * pts[i].z;

        mxx += nx * nx;
        mxy += nx * ny;
        mxz += nx * nz;
        myy += ny * ny;
        myz += ny * nz;
        mzz += nz * nz;

        yx += dot * nx;
        yy += dot * ny;
        yz += dot * nz;
    }

    // Add regularization parameter to the diagonal (Ridge Regression / Tikhonov)
    double lambda = regularization;
    mxx += lambda;
    myy += lambda;
    mzz += lambda;

    yx += lambda * fallback.x;
    yy += lambda * fallback.y;
    yz += lambda * fallback.z;

    double det = mxx * (myy * mzz - myz * myz) - mxy * (mxy * mzz - myz * mxz) + mxz * (mxy * myz - myy * mxz);

    if (Math::abs(det) < 1e-7) {
        return fallback;
    }

    double idet = 1.0 / det;
    double rxx = (myy * mzz - myz * myz) * idet;
    double rxy = (myz * mxz - mxy * mzz) * idet;
    double rxz = (mxy * myz - myy * mxz) * idet;
    double ryy = (mxx * mzz - mxz * mxz) * idet;
    double ryz = (mxy * mxz - mxx * myz) * idet;
    double rzz = (mxx * myy - mxy * mxy) * idet;

    Vector3 sol;
    sol.x = rxx * yx + rxy * yy + rxz * yz;
    sol.y = rxy * yx + ryy * yy + ryz * yz;
    sol.z = rxz * yx + ryz * yy + rzz * yz;

    // Clamp to cell bounds to keep solution inside voxel
    if (sol.x < cell_min.x || sol.x > cell_max.x ||
        sol.y < cell_min.y || sol.y > cell_max.y ||
        sol.z < cell_min.z || sol.z > cell_max.z) {
        sol.x = Math::clamp(sol.x, (float)cell_min.x, (float)cell_max.x);
        sol.y = Math::clamp(sol.y, (float)cell_min.y, (float)cell_max.y);
        sol.z = Math::clamp(sol.z, (float)cell_min.z, (float)cell_max.z);
    }

    return sol;
}

void DCChunk::generate_mesh_from_density_grid() {
    _clear_collision();
    _clear_occluder();
    _clear_liquid();

    if (!density_grid.is_valid()) {
        Ref<ArrayMesh> current_mesh = get_mesh();
        if (current_mesh.is_valid()) current_mesh->clear_surfaces();
        return;
    }

    float threshold = density_grid->get_surface_threshold();
    int grid_dim_x = density_grid->get_grid_size_x();
    int grid_dim_y = density_grid->get_grid_size_y();
    int grid_dim_z = density_grid->get_grid_size_z();
    const PackedFloat32Array &density_array = density_grid->get_density_data();
    const PackedByteArray &material_array = density_grid->get_material_data();
    const float *density_data = density_array.ptr();
    const uint8_t *material_data = material_array.ptr();

    auto sample_density = [&](int wx, int wy, int wz) -> float {
        if (wx < 0 || wx >= grid_dim_x || wy < 0 || wy >= grid_dim_y || wz < 0 || wz >= grid_dim_z) {
            if (wy >= grid_dim_y) {
                return -1.0f; // Sky above the grid is empty/air
            }
            if (wy < 0) {
                return 1.0f; // Bedrock below the grid is solid
            }
            int clamped_x = wx < 0 ? 0 : (wx >= grid_dim_x ? grid_dim_x - 1 : wx);
            int clamped_z = wz < 0 ? 0 : (wz >= grid_dim_z ? grid_dim_z - 1 : wz);
            int idx = clamped_x + grid_dim_x * (wy + grid_dim_y * clamped_z);
            return density_data[idx];
        }
        int idx = wx + grid_dim_x * (wy + grid_dim_y * wz);
        return density_data[idx];
    };

    auto sample_material = [&](int wx, int wy, int wz) -> int {
        if (wx < 0 || wx >= grid_dim_x || wy < 0 || wy >= grid_dim_y || wz < 0 || wz >= grid_dim_z) {
            return 0;
        }
        int idx = wx + grid_dim_x * (wy + grid_dim_y * wz);
        int mat = (int)material_data[idx];
        if (mat == liquid_material_id || (mat >= 255 - flow_spread_limit && mat <= 255)) {
            return 0;
        }
        return mat;
    };

    auto get_grid_gradient = [&](const Vector3i &pos) -> Vector3 {
        float dx = sample_density(pos.x + 1, pos.y, pos.z) - sample_density(pos.x - 1, pos.y, pos.z);
        float dy = sample_density(pos.x, pos.y + 1, pos.z) - sample_density(pos.x, pos.y - 1, pos.z);
        float dz = sample_density(pos.x, pos.y, pos.z + 1) - sample_density(pos.x, pos.y, pos.z - 1);
        Vector3 g(dx, dy, dz);
        if (g.length_squared() < 1e-8f) {
            return Vector3(0, 0, 0);
        }
        return -g.normalized();
    };

    const Vector3i corner_offsets[8] = {
        Vector3i(0, 0, 0), Vector3i(1, 0, 0), Vector3i(1, 1, 0), Vector3i(0, 1, 0),
        Vector3i(0, 0, 1), Vector3i(1, 0, 1), Vector3i(1, 1, 1), Vector3i(0, 1, 1)
    };

    const int cube_edges[12][2] = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0},
        {4, 5}, {5, 6}, {6, 7}, {7, 4},
        {0, 4}, {1, 5}, {2, 6}, {3, 7}
    };

    std::unordered_map<Vector3i, CellVertexInfo, Vector3iHash> cell_vertices;

    // Pass 1: Find all active cells in the chunk region (including boundary stitching cells)
    for (int cz = -1; cz < chunk_size; ++cz) {
        for (int cy = -1; cy < chunk_size; ++cy) {
            for (int cx = -1; cx < chunk_size; ++cx) {
                Vector3i local_cell(cx, cy, cz);
                Vector3i world_cell = chunk_grid_offset + local_cell;

                float values[8];
                bool is_solid[8];
                for (int i = 0; i < 8; ++i) {
                    Vector3i cp = world_cell + corner_offsets[i];
                    values[i] = sample_density(cp.x, cp.y, cp.z);
                    is_solid[i] = values[i] > threshold;
                }

                std::vector<Vector3> pts;
                std::vector<Vector3> norms;
                pts.reserve(12);
                norms.reserve(12);

                for (int j = 0; j < 12; ++j) {
                    int c1 = cube_edges[j][0];
                    int c2 = cube_edges[j][1];
                    if (is_solid[c1] != is_solid[c2]) {
                        Vector3i cp1 = world_cell + corner_offsets[c1];
                        Vector3i cp2 = world_cell + corner_offsets[c2];

                        float val1 = values[c1];
                        float val2 = values[c2];

                        float t = 0.5f;
                        if (Math::abs(val1 - val2) > 1e-6f) {
                            t = (threshold - val1) / (val2 - val1);
                        }

                        Vector3 p_grid = Vector3(cp1) + t * Vector3(cp2 - cp1);
                        pts.push_back(p_grid);

                        Vector3 n1 = get_grid_gradient(cp1);
                        Vector3 n2 = get_grid_gradient(cp2);
                        Vector3 n_interp = n1 + t * (n2 - n1);
                        if (n_interp.length_squared() > 1e-8f) {
                            norms.push_back(n_interp.normalized());
                        } else {
                            norms.push_back(Vector3(0, 1, 0));
                        }
                    }
                }

                if (pts.empty()) continue; // Inactive cell

                Vector3 cell_min = Vector3(world_cell);
                Vector3 cell_max = Vector3(world_cell + Vector3i(1, 1, 1));

                Vector3 vertex_grid;
                if (use_qef) {
                    vertex_grid = solve_qef_cramer(pts, norms, cell_min, cell_max, qef_regularization);
                } else {
                    Vector3 sum(0, 0, 0);
                    for (const Vector3 &p : pts) sum += p;
                    vertex_grid = sum / (float)pts.size();
                }

                Vector3 vertex_local = (vertex_grid - Vector3(chunk_grid_offset)) * voxel_size;

                Vector3 normal_avg(0, 0, 0);
                for (const Vector3 &n : norms) normal_avg += n;
                if (normal_avg.length_squared() > 1e-8f) {
                    normal_avg = normal_avg.normalized();
                } else {
                    normal_avg = Vector3(0, 1, 0);
                }

                int mat_idx = -1;
                float max_solid_density = -99999.0f;
                for (int i = 0; i < 8; ++i) {
                    if (is_solid[i]) {
                        Vector3i cp = world_cell + corner_offsets[i];
                        int mat = sample_material(cp.x, cp.y, cp.z);
                        if (mat > 0) {
                            if (values[i] > max_solid_density) {
                                max_solid_density = values[i];
                                mat_idx = mat;
                            }
                        }
                    }
                }

                if (mat_idx == -1) {
                    int min_density_idx = 0;
                    float min_density = values[0];
                    for (int i = 1; i < 8; ++i) {
                        if (values[i] < min_density) {
                            min_density = values[i];
                            min_density_idx = i;
                        }
                    }
                    Vector3i cp = world_cell + corner_offsets[min_density_idx];
                    mat_idx = sample_material(cp.x, cp.y, cp.z);
                }

                CellVertexInfo info;
                info.position = vertex_local;
                info.normal = normal_avg;
                info.material_id = mat_idx;
                cell_vertices[world_cell] = info;
            }
        }
    }

    PackedVector3Array master_vertices;
    PackedVector3Array master_normals;
    std::vector<int> master_triangle_materials;
    PackedInt32Array master_indices;

    auto get_or_create_vertex = [&](const Vector3i &world_cell) -> int {
        auto it = cell_vertices.find(world_cell);
        if (it == cell_vertices.end()) return -1;

        CellVertexInfo &info = it->second;
        if (info.master_index != -1) return info.master_index;

        int idx = master_vertices.size();
        master_vertices.append(info.position);
        master_normals.append(info.normal);
        info.master_index = idx;
        return idx;
    };

    // Pass 2: Quad generation for edges inside the chunk bounds
    for (int z = 0; z < chunk_size; ++z) {
        for (int y = 0; y < chunk_size; ++y) {
            for (int x = 0; x < chunk_size; ++x) {
                Vector3i local_pos(x, y, z);
                Vector3i world_pos = chunk_grid_offset + local_pos;

                // X edge
                if (world_pos.x + 1 < grid_dim_x) {
                    float d_start = sample_density(world_pos.x, world_pos.y, world_pos.z);
                    float d_end = sample_density(world_pos.x + 1, world_pos.y, world_pos.z);
                    bool start_solid = d_start > threshold;
                    bool end_solid = d_end > threshold;

                    if (start_solid != end_solid) {
                        Vector3i c0 = world_pos;
                        Vector3i c1 = world_pos - Vector3i(0, 0, 1);
                        Vector3i c2 = world_pos - Vector3i(0, 1, 1);
                        Vector3i c3 = world_pos - Vector3i(0, 1, 0);

                        int i0 = get_or_create_vertex(c0);
                        int i1 = get_or_create_vertex(c1);
                        int i2 = get_or_create_vertex(c2);
                        int i3 = get_or_create_vertex(c3);

                        if (i0 != -1 && i1 != -1 && i2 != -1 && i3 != -1) {
                            int mat = start_solid ? sample_material(world_pos.x, world_pos.y, world_pos.z)
                                                  : sample_material(world_pos.x + 1, world_pos.y, world_pos.z);

                            if (start_solid) {
                                master_indices.append(i0); master_indices.append(i1); master_indices.append(i2);
                                master_indices.append(i0); master_indices.append(i2); master_indices.append(i3);
                            } else {
                                master_indices.append(i0); master_indices.append(i3); master_indices.append(i2);
                                master_indices.append(i0); master_indices.append(i2); master_indices.append(i1);
                            }
                            master_triangle_materials.push_back(mat);
                            master_triangle_materials.push_back(mat);
                        }
                    }
                }

                // Y edge
                if (world_pos.y + 1 < grid_dim_y) {
                    float d_start = sample_density(world_pos.x, world_pos.y, world_pos.z);
                    float d_end = sample_density(world_pos.x, world_pos.y + 1, world_pos.z);
                    bool start_solid = d_start > threshold;
                    bool end_solid = d_end > threshold;

                    if (start_solid != end_solid) {
                        Vector3i c0 = world_pos;
                        Vector3i c1 = world_pos - Vector3i(1, 0, 0);
                        Vector3i c2 = world_pos - Vector3i(1, 0, 1);
                        Vector3i c3 = world_pos - Vector3i(0, 0, 1);

                        int i0 = get_or_create_vertex(c0);
                        int i1 = get_or_create_vertex(c1);
                        int i2 = get_or_create_vertex(c2);
                        int i3 = get_or_create_vertex(c3);

                        if (i0 != -1 && i1 != -1 && i2 != -1 && i3 != -1) {
                            int mat = start_solid ? sample_material(world_pos.x, world_pos.y, world_pos.z)
                                                  : sample_material(world_pos.x, world_pos.y + 1, world_pos.z);

                            if (start_solid) {
                                master_indices.append(i0); master_indices.append(i1); master_indices.append(i2);
                                master_indices.append(i0); master_indices.append(i2); master_indices.append(i3);
                            } else {
                                master_indices.append(i0); master_indices.append(i3); master_indices.append(i2);
                                master_indices.append(i0); master_indices.append(i2); master_indices.append(i1);
                            }
                            master_triangle_materials.push_back(mat);
                            master_triangle_materials.push_back(mat);
                        }
                    }
                }

                // Z edge
                if (world_pos.z + 1 < grid_dim_z) {
                    float d_start = sample_density(world_pos.x, world_pos.y, world_pos.z);
                    float d_end = sample_density(world_pos.x, world_pos.y, world_pos.z + 1);
                    bool start_solid = d_start > threshold;
                    bool end_solid = d_end > threshold;

                    if (start_solid != end_solid) {
                        Vector3i c0 = world_pos;
                        Vector3i c1 = world_pos - Vector3i(0, 1, 0);
                        Vector3i c2 = world_pos - Vector3i(1, 1, 0);
                        Vector3i c3 = world_pos - Vector3i(1, 0, 0);

                        int i0 = get_or_create_vertex(c0);
                        int i1 = get_or_create_vertex(c1);
                        int i2 = get_or_create_vertex(c2);
                        int i3 = get_or_create_vertex(c3);

                        if (i0 != -1 && i1 != -1 && i2 != -1 && i3 != -1) {
                            int mat = start_solid ? sample_material(world_pos.x, world_pos.y, world_pos.z)
                                                  : sample_material(world_pos.x, world_pos.y, world_pos.z + 1);

                            if (start_solid) {
                                master_indices.append(i0); master_indices.append(i1); master_indices.append(i2);
                                master_indices.append(i0); master_indices.append(i2); master_indices.append(i3);
                            } else {
                                master_indices.append(i0); master_indices.append(i3); master_indices.append(i2);
                                master_indices.append(i0); master_indices.append(i2); master_indices.append(i1);
                            }
                            master_triangle_materials.push_back(mat);
                            master_triangle_materials.push_back(mat);
                        }
                    }
                }
            }
        }
    }

    if (master_vertices.is_empty()) return;

    if (flip_normals) {
        for (int i = 0; i < master_normals.size(); ++i) {
            master_normals[i] = -master_normals[i];
        }
    }

    // Split into Material Surfaces
    std::map<int, SurfaceBuilder> surfaces;
    int num_triangles = master_triangle_materials.size();
    for (int t = 0; t < num_triangles; ++t) {
        int mat_id = master_triangle_materials[t];
        SurfaceBuilder &surf = surfaces[mat_id];

        for (int k = 0; k < 3; ++k) {
            int master_idx = master_indices[t * 3 + k];
            auto it = surf.index_cache.find(master_idx);
            if (it != surf.index_cache.end()) {
                surf.indices.append(it->second);
            } else {
                int new_surf_idx = surf.vertices.size();
                surf.vertices.append(master_vertices[master_idx]);
                surf.normals.append(master_normals[master_idx]);
                surf.indices.append(new_surf_idx);
                surf.index_cache[master_idx] = new_surf_idx;
            }
        }
    }

    Ref<ArrayMesh> array_mesh = get_mesh();
    if (!array_mesh.is_valid()) {
        array_mesh.instantiate();
        set_mesh(array_mesh);
    }
    array_mesh->clear_surfaces();

    PackedVector3Array total_collision_vertices;
    PackedInt32Array total_collision_indices;
    int col_offset = 0;

    for (auto const& [mat_id, builder] : surfaces) {
        if (builder.vertices.is_empty()) continue;

        Array arrays;
        arrays.resize(Mesh::ARRAY_MAX);
        arrays[Mesh::ARRAY_VERTEX] = builder.vertices;
        arrays[Mesh::ARRAY_NORMAL] = builder.normals;
        arrays[Mesh::ARRAY_INDEX] = builder.indices;

        PackedVector2Array uvs;
        uvs.resize(builder.vertices.size());
        float uv_scale = voxel_size > 0 ? voxel_size : 1.0f;
        for (int i = 0; i < builder.vertices.size(); ++i) {
            uvs[i] = Vector2(builder.vertices[i].x / uv_scale, builder.vertices[i].y / uv_scale);
        }
        arrays[Mesh::ARRAY_TEX_UV] = uvs;

        array_mesh->add_surface_from_arrays(Mesh::PRIMITIVE_TRIANGLES, arrays);

        if (mat_id >= 0 && mat_id < materials.size()) {
            Ref<Material> mat = materials[mat_id];
            if (mat.is_valid()) {
                set_surface_override_material(array_mesh->get_surface_count() - 1, mat);
            }
        }

        if (generate_collision || generate_occluder) {
            total_collision_vertices.append_array(builder.vertices);
            for (int idx : builder.indices) {
                total_collision_indices.append(idx + col_offset);
            }
            col_offset += builder.vertices.size();
        }
    }

    if (generate_collision) _generate_collision(total_collision_vertices, total_collision_indices);
    if (generate_occluder) _generate_occluder(total_collision_vertices, total_collision_indices);
    if (liquid_material.is_valid()) {
        _generate_liquid_mesh();
    }
}

void DCChunk::_generate_collision(const PackedVector3Array &p_vertices, const PackedInt32Array &p_indices) {
    if (p_vertices.is_empty() || p_indices.is_empty()) {
        _clear_collision();
        return;
    }

    StaticBody3D *static_body = memnew(StaticBody3D);
    add_child(static_body);
    static_body->set_owner(get_owner());

    collision_shape.instantiate();
    if (!collision_shape.is_valid()) {
        static_body->queue_free();
        return;
    }

    PackedVector3Array collision_face_vertices;
    collision_face_vertices.resize(p_indices.size());

    for (int i = 0; i < p_indices.size(); ++i) {
        int vertex_index = p_indices[i];
        if (vertex_index >= 0 && vertex_index < p_vertices.size()) {
            collision_face_vertices[i] = p_vertices[vertex_index];
        } else {
            collision_face_vertices.clear();
            static_body->queue_free();
            collision_shape.unref();
            return;
        }
    }
    collision_shape->set_faces(collision_face_vertices);

    CollisionShape3D *cs = memnew(CollisionShape3D);
    static_body->add_child(cs);
    cs->set_owner(get_owner());
    cs->set_shape(collision_shape);
}

void DCChunk::_clear_collision() {
    for (int i = get_child_count() - 1; i >= 0; --i) {
        Node *child = get_child(i);
        if (Object::cast_to<StaticBody3D>(child)) {
            child->queue_free();
        }
    }
    collision_shape.unref();
}

void DCChunk::_generate_occluder(const PackedVector3Array &p_vertices, const PackedInt32Array &p_indices) {
    if (p_vertices.is_empty() || p_indices.is_empty()) {
        _clear_occluder();
        return;
    }

    OccluderInstance3D *occluder_instance = memnew(OccluderInstance3D);
    add_child(occluder_instance);
    occluder_instance->set_owner(get_owner());

    occluder_shape.instantiate();
    if (!occluder_shape.is_valid()) {
        occluder_instance->queue_free();
        return;
    }

    occluder_shape->set_vertices(p_vertices);
    occluder_shape->set_indices(p_indices);
    occluder_instance->set_occluder(occluder_shape);
}

void DCChunk::_clear_occluder() {
    for (int i = get_child_count() - 1; i >= 0; --i) {
        Node *child = get_child(i);
        if (Object::cast_to<OccluderInstance3D>(child)) {
            child->queue_free();
        }
    }
    occluder_shape.unref();
}

void DCChunk::_generate_liquid_mesh() {
    if (!density_grid.is_valid()) return;

    int dim_x = density_grid->get_grid_size_x();
    int dim_y = density_grid->get_grid_size_y();
    int dim_z = density_grid->get_grid_size_z();
    float surf_thresh = density_grid->get_surface_threshold();

    bool has_liquid = false;
    for (int z = 0; z < chunk_size; ++z) {
        for (int y = 0; y < chunk_size; ++y) {
            for (int x = 0; x < chunk_size; ++x) {
                Vector3i world_pos = chunk_grid_offset + Vector3i(x, y, z);
                if (world_pos.x >= 0 && world_pos.x < dim_x &&
                    world_pos.y >= 0 && world_pos.y < dim_y &&
                    world_pos.z >= 0 && world_pos.z < dim_z) {
                    int mat = density_grid->get_material_id(world_pos);
                    if (mat == liquid_material_id || (mat >= 255 - flow_spread_limit && mat <= 255)) {
                        has_liquid = true;
                        break;
                    }
                }
            }
            if (has_liquid) break;
        }
        if (has_liquid) break;
    }

    if (!has_liquid) return;

    auto get_liquid_density = [&](const Vector3i &wpos) -> float {
        if (wpos.x < 0 || wpos.x >= dim_x || wpos.y < 0 || wpos.y >= dim_y || wpos.z < 0 || wpos.z >= dim_z) {
            return -1.0f;
        }
        if (density_grid->get_cell(wpos) > surf_thresh) {
            return -1.0f;
        }
        int mat = density_grid->get_material_id(wpos);
        if (mat == liquid_material_id || (mat >= 255 - flow_spread_limit && mat <= 255)) {
            int max_level = flow_spread_limit + 1;
            int level = max_level;
            if (mat >= 255 - flow_spread_limit && mat <= 255) {
                level = mat - (255 - max_level);
            }
            float density = 0.05f + 0.95f * (float)(level - 1) / (float)flow_spread_limit;
            return density;
        }
        return -1.0f;
    };

    PackedVector3Array output_vertices;
    PackedInt32Array output_triangles;
    std::unordered_map<uint64_t, int> liquid_vertex_cache;

    const Vector3i c_offsets[8] = {
        Vector3i(0, 0, 0), Vector3i(1, 0, 0), Vector3i(1, 1, 0), Vector3i(0, 1, 0),
        Vector3i(0, 0, 1), Vector3i(1, 0, 1), Vector3i(1, 1, 1), Vector3i(0, 1, 1)
    };

    for (int z = 0; z < chunk_size; ++z) {
        for (int y = 0; y < chunk_size; ++y) {
            for (int x = 0; x < chunk_size; ++x) {
                Vector3i local_pos(x, y, z);
                Vector3i world_pos_base = chunk_grid_offset + local_pos;

                float corner_values[8];
                for (int i = 0; i < 8; ++i) {
                    corner_values[i] = get_liquid_density(world_pos_base + c_offsets[i]);
                }

                int cube_index = 0;
                for (int i = 0; i < 8; ++i) {
                    if (corner_values[i] > 0.0f) {
                        cube_index |= (1 << i);
                    }
                }

                const int* tri_table_row = McTables::TRI_TABLE[cube_index];
                if (tri_table_row[0] == -1) continue;

                Vector3 corner_locations_local[8];
                for (int i = 0; i < 8; ++i) {
                    corner_locations_local[i] = Vector3(local_pos + c_offsets[i]) * voxel_size;
                }

                int triangle_indices[3];
                for (int i = 0; tri_table_row[i] != -1; i += 3) {
                    for (int j = 0; j < 3; ++j) {
                        int edge_index = tri_table_row[i + j];
                        int corner_a_idx = McTables::CORNER_INDEX_A_FROM_EDGE[edge_index];
                        int corner_b_idx = McTables::CORNER_INDEX_B_FROM_EDGE[edge_index];
                        
                        int cached_idx = -1;
                        uint64_t edge_key = 0;
                        if (smooth_normals) {
                            Vector3i world_corner_a = world_pos_base + c_offsets[corner_a_idx];
                            Vector3i world_corner_b = world_pos_base + c_offsets[corner_b_idx];
                            int px = world_corner_a.x < world_corner_b.x ? world_corner_a.x : world_corner_b.x;
                            int py = world_corner_a.y < world_corner_b.y ? world_corner_a.y : world_corner_b.y;
                            int pz = world_corner_a.z < world_corner_b.z ? world_corner_a.z : world_corner_b.z;
                            int axis = 0;
                            if (world_corner_a.y != world_corner_b.y) axis = 1;
                            else if (world_corner_a.z != world_corner_b.z) axis = 2;

                            edge_key = (((uint64_t)px & 0xFFFFF) << 42) |
                                       (((uint64_t)py & 0xFFFFF) << 22) |
                                       (((uint64_t)pz & 0xFFFFF) << 2) |
                                       ((uint64_t)axis & 0x3);

                            auto cache_it = liquid_vertex_cache.find(edge_key);
                            if (cache_it != liquid_vertex_cache.end()) {
                                cached_idx = cache_it->second;
                            }
                        }

                        if (cached_idx != -1) {
                            triangle_indices[j] = cached_idx;
                        } else {
                            float val1 = corner_values[corner_a_idx];
                            float val2 = corner_values[corner_b_idx];
                            Vector3 p1 = corner_locations_local[corner_a_idx];
                            Vector3 p2 = corner_locations_local[corner_b_idx];
                            Vector3 vert_pos;
                            if (Math::abs(val1 - val2) < CMP_EPSILON) {
                                vert_pos = (p1 + p2) * 0.5f;
                            } else {
                                float t = (0.0f - val1) / (val2 - val1);
                                vert_pos = p1 + t * (p2 - p1);
                            }

                            int new_vertex_index = output_vertices.size();
                            output_vertices.append(vert_pos);
                            if (smooth_normals) {
                                liquid_vertex_cache[edge_key] = new_vertex_index;
                            }
                            triangle_indices[j] = new_vertex_index;
                        }
                    }

                    output_triangles.append(triangle_indices[0]);
                    output_triangles.append(triangle_indices[1]);
                    output_triangles.append(triangle_indices[2]);
                }
            }
        }
    }

    if (output_vertices.is_empty()) return;

    PackedVector3Array output_normals;
    output_normals.resize(output_vertices.size());
    for(int i = 0; i < output_normals.size(); ++i) {
        output_normals[i] = Vector3(0, 0, 0);
    }
    for (int i = 0; i < output_triangles.size(); i += 3) {
        int i1 = output_triangles[i];
        int i2 = output_triangles[i+1];
        int i3 = output_triangles[i+2];

        Vector3 v1 = output_vertices[i1];
        Vector3 v2 = output_vertices[i2];
        Vector3 v3 = output_vertices[i3];

        Vector3 face_normal = (v2 - v1).cross(v3 - v1);
        output_normals[i1] += face_normal;
        output_normals[i2] += face_normal;
        output_normals[i3] += face_normal;
    }
    float normal_sign = flip_normals ? -1.0f : 1.0f;
    for(int i = 0; i < output_normals.size(); ++i) {
        if (output_normals[i].length_squared() < 1e-8f) {
            output_normals[i] = Vector3(0, 1, 0) * normal_sign;
        } else {
            output_normals[i] = output_normals[i].normalized() * normal_sign;
        }
    }

    PackedVector2Array output_uvs;
    output_uvs.resize(output_vertices.size());
    float uv_scale = voxel_size > 0 ? voxel_size : 1.0f;
    for(int i = 0; i < output_vertices.size(); ++i) {
        output_uvs[i] = Vector2(output_vertices[i].x / uv_scale, output_vertices[i].z / uv_scale);
    }

    Ref<ArrayMesh> liquid_arr_mesh;
    liquid_arr_mesh.instantiate();

    Array arrays;
    arrays.resize(Mesh::ARRAY_MAX);
    arrays[Mesh::ARRAY_VERTEX] = output_vertices;
    arrays[Mesh::ARRAY_NORMAL] = output_normals;
    arrays[Mesh::ARRAY_INDEX] = output_triangles;
    arrays[Mesh::ARRAY_TEX_UV] = output_uvs;

    liquid_arr_mesh->add_surface_from_arrays(Mesh::PRIMITIVE_TRIANGLES, arrays);

    MeshInstance3D* liquid_mesh_node = memnew(MeshInstance3D);
    liquid_mesh_node->set_name("LiquidMesh");
    add_child(liquid_mesh_node);
    liquid_mesh_node->set_owner(get_owner());
    liquid_mesh_node->set_mesh(liquid_arr_mesh);
    if (liquid_material.is_valid()) {
        liquid_mesh_node->set_material_override(liquid_material);
    }

    if (generate_liquid_trigger) {
        Area3D* area = memnew(Area3D);
        area->set_name("LiquidTrigger");
        add_child(area);
        area->set_owner(get_owner());

        CollisionShape3D* shape_node = memnew(CollisionShape3D);
        shape_node->set_name("CollisionShape3D");
        area->add_child(shape_node);
        shape_node->set_owner(get_owner());

        Ref<ConcavePolygonShape3D> liquid_trigger_shape;
        liquid_trigger_shape.instantiate();

        PackedVector3Array faces;
        faces.resize(output_triangles.size());
        for (int i = 0; i < output_triangles.size(); ++i) {
            faces[i] = output_vertices[output_triangles[i]];
        }
        liquid_trigger_shape->set_faces(faces);
        shape_node->set_shape(liquid_trigger_shape);
    }
}

void DCChunk::_clear_liquid() {
    for (int i = get_child_count() - 1; i >= 0; --i) {
        Node *child = get_child(i);
        if (child && (child->get_name() == String("LiquidMesh") || child->get_name() == String("LiquidTrigger"))) {
            child->queue_free();
        }
    }
}

void DCChunk::set_chunk_size(int p_size) { chunk_size = p_size > 0 ? p_size : 1; }
int DCChunk::get_chunk_size() const { return chunk_size; }

void DCChunk::set_voxel_size(float p_size) { voxel_size = p_size > 0 ? p_size : 0.001f; }
float DCChunk::get_voxel_size() const { return voxel_size; }

void DCChunk::set_chunk_grid_offset(const Vector3i &p_offset) { chunk_grid_offset = p_offset; }
Vector3i DCChunk::get_chunk_grid_offset() const { return chunk_grid_offset; }

void DCChunk::set_density_grid(const Ref<DensityGrid> &p_grid) { density_grid = p_grid; }
Ref<DensityGrid> DCChunk::get_density_grid() const { return density_grid; }

void DCChunk::set_generate_collision(bool p_generate) { generate_collision = p_generate; }
bool DCChunk::get_generate_collision() const { return generate_collision; }

void DCChunk::set_generate_occluder(bool p_generate) { generate_occluder = p_generate; }
bool DCChunk::get_generate_occluder() const { return generate_occluder; }

void DCChunk::set_materials(const TypedArray<Material> &p_materials) { materials = p_materials; }
TypedArray<Material> DCChunk::get_materials() const { return materials; }

void DCChunk::set_compute_shader(const Ref<RDShaderFile> &p_shader) { compute_shader = p_shader; }
Ref<RDShaderFile> DCChunk::get_compute_shader() const { return compute_shader; }

void DCChunk::set_liquid_material(const Ref<Material> &p_material) { liquid_material = p_material; }
Ref<Material> DCChunk::get_liquid_material() const { return liquid_material; }

void DCChunk::set_liquid_material_id(int p_id) { liquid_material_id = Math::clamp(p_id, 0, 255); }
int DCChunk::get_liquid_material_id() const { return liquid_material_id; }

void DCChunk::set_generate_liquid_trigger(bool p_enabled) { generate_liquid_trigger = p_enabled; }
bool DCChunk::get_generate_liquid_trigger() const { return generate_liquid_trigger; }

void DCChunk::set_flow_spread_limit(int p_limit) { flow_spread_limit = p_limit > 0 ? p_limit : 1; }
int DCChunk::get_flow_spread_limit() const { return flow_spread_limit; }

void DCChunk::set_smooth_normals(bool p_smooth) { smooth_normals = p_smooth; }
bool DCChunk::get_smooth_normals() const { return smooth_normals; }

void DCChunk::set_flip_normals(bool p_flip) { flip_normals = p_flip; }
bool DCChunk::get_flip_normals() const { return flip_normals; }

void DCChunk::set_use_qef(bool p_use) { use_qef = p_use; }
bool DCChunk::get_use_qef() const { return use_qef; }

void DCChunk::set_qef_regularization(float p_reg) { qef_regularization = p_reg; }
float DCChunk::get_qef_regularization() const { return qef_regularization; }

} // namespace godot
