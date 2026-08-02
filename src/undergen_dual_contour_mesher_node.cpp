// src/undergen_dual_contour_mesher_node.cpp
#include "undergen_dual_contour_mesher_node.h"
#include "dc_chunk.h"
#include <godot_cpp/variant/utility_functions.hpp>
#include <godot_cpp/core/math.hpp>
#include <godot_cpp/classes/geometry_instance3d.hpp>

namespace godot {

UnderGenDualContourMesherNode::UnderGenDualContourMesherNode() {}
UnderGenDualContourMesherNode::~UnderGenDualContourMesherNode() {}

void UnderGenDualContourMesherNode::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_use_qef", "use"), &UnderGenDualContourMesherNode::set_use_qef);
    ClassDB::bind_method(D_METHOD("get_use_qef"), &UnderGenDualContourMesherNode::get_use_qef);
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "use_qef"), "set_use_qef", "get_use_qef");

    ClassDB::bind_method(D_METHOD("set_qef_regularization", "regularization"), &UnderGenDualContourMesherNode::set_qef_regularization);
    ClassDB::bind_method(D_METHOD("get_qef_regularization"), &UnderGenDualContourMesherNode::get_qef_regularization);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "qef_regularization", PROPERTY_HINT_RANGE, "0.0, 1.0, 0.00001, exp, or_greater"), "set_qef_regularization", "get_qef_regularization");

    ClassDB::bind_method(D_METHOD("set_stepped_transitions", "stepped"), &UnderGenDualContourMesherNode::set_stepped_transitions);
    ClassDB::bind_method(D_METHOD("get_stepped_transitions"), &UnderGenDualContourMesherNode::get_stepped_transitions);
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "stepped_transitions"), "set_stepped_transitions", "get_stepped_transitions");
}

void UnderGenDualContourMesherNode::set_use_qef(bool p_use) { use_qef = p_use; }
bool UnderGenDualContourMesherNode::get_use_qef() const { return use_qef; }

void UnderGenDualContourMesherNode::set_qef_regularization(float p_reg) { qef_regularization = p_reg; }
float UnderGenDualContourMesherNode::get_qef_regularization() const { return qef_regularization; }

void UnderGenDualContourMesherNode::set_stepped_transitions(bool p_stepped) { stepped_transitions = p_stepped; }
bool UnderGenDualContourMesherNode::get_stepped_transitions() const { return stepped_transitions; }

void UnderGenDualContourMesherNode::execute_with_parent(const Dictionary &inputs, Dictionary &outputs, Node3D* parent_node) {
    if (!parent_node) {
        UtilityFunctions::printerr("UnderGenDualContourMesherNode: No parent node provided.");
        return;
    }

    Dictionary context = inputs.get(0, Dictionary());
    if (context.is_empty()) return;

    Ref<DensityGrid> grid = context.get("grid", Ref<DensityGrid>());
    if (grid.is_null()) {
        UtilityFunctions::printerr("UnderGenDualContourMesherNode: Grid not found in context.");
        return;
    }

    Vector3i dims = grid->get_grid_dimensions();
    int chunk_sz = get_chunk_size();
    float vox_sz = get_voxel_size();
    bool gen_col = get_generate_collision();
    bool gen_occ = get_generate_occluder();
    TypedArray<Material> terrain_mats = get_terrain_materials();
    Ref<Material> liquid_mat = get_liquid_material();
    int liquid_mat_id = get_liquid_material_id();
    bool gen_liquid_trig = get_generate_liquid_trigger();
    Ref<RDShaderFile> comp_shader = get_compute_shader();
    bool smooth_norms = get_smooth_normals();
    bool flip_norms = get_flip_normals();
    bool cast_shad = get_cast_shadows();

    int count_x = Math::max(1, (int)Math::ceil((float)(dims.x - 1) / chunk_sz));
    int count_y = Math::max(1, (int)Math::ceil((float)(dims.y - 1) / chunk_sz));
    int count_z = Math::max(1, (int)Math::ceil((float)(dims.z - 1) / chunk_sz));

    UtilityFunctions::print("UnderGenDualContourMesherNode: Spawning ", count_x * count_y * count_z, " dual contouring chunks...");

    for (int x = 0; x < count_x; ++x) {
        for (int y = 0; y < count_y; ++y) {
            for (int z = 0; z < count_z; ++z) {
                DCChunk* chunk = memnew(DCChunk);
                chunk->set_name("TerrainChunk_DC_" + String::num_int64(x) + "_" + String::num_int64(y) + "_" + String::num_int64(z));
                chunk->set_chunk_size(chunk_sz);
                chunk->set_voxel_size(vox_sz);
                chunk->set_chunk_grid_offset(Vector3i(x, y, z) * chunk_sz);
                chunk->set_density_grid(grid);
                chunk->set_generate_collision(gen_col);
                chunk->set_generate_occluder(gen_occ);
                chunk->set_smooth_normals(smooth_norms);
                chunk->set_flip_normals(flip_norms);
                chunk->set_cast_shadows_setting(cast_shad ? GeometryInstance3D::SHADOW_CASTING_SETTING_ON : GeometryInstance3D::SHADOW_CASTING_SETTING_OFF);
                chunk->set_use_qef(use_qef);
                chunk->set_qef_regularization(qef_regularization);
                chunk->call("set_stepped_transitions", stepped_transitions);

                if (liquid_mat.is_valid()) {
                    chunk->call("set_liquid_material", liquid_mat);
                    chunk->call("set_liquid_material_id", liquid_mat_id);
                    chunk->call("set_generate_liquid_trigger", gen_liquid_trig);
                    
                    int flow_spread_limit = context.get("flow_spread_limit", 7);
                    chunk->call("set_flow_spread_limit", flow_spread_limit);
                }
                if (!terrain_mats.is_empty()) {
                    chunk->set_materials(terrain_mats);
                }
                if (context.has("material_thicknesses")) {
                    chunk->call("set_material_thicknesses", context["material_thicknesses"]);
                }
                if (context.has("material_stepped")) {
                    chunk->call("set_material_stepped", context["material_stepped"]);
                }
                if (comp_shader.is_valid()) {
                    chunk->set_compute_shader(comp_shader);
                }
                chunk->set_position(Vector3(x, y, z) * chunk_sz * vox_sz);
                parent_node->add_child(chunk);
                chunk->call_deferred("generate_mesh_from_density_grid");
            }
        }
    }

    outputs[0] = context;
    UtilityFunctions::print("UnderGenDualContourMesherNode: Chunk spawning complete.");
}

} // namespace godot
