// src/register_types.cpp
#include "register_types.h"

// Existing classes
#include "density_grid.h"
#include "mc_chunk.h"

// New pipeline infrastructure
#include "undergen_node.h"
#include "undergen_point_set.h"
#include "undergen_pipeline.h"
#include "undergen_world_3d.h"
#include "level_grammar_spec.h"

// Concrete pipeline nodes
#include "undergen_bsp_placer_node.h"
#include "undergen_astar_carver_node.h"
#include "undergen_bezier_carver_node.h"
#include "undergen_noise_node.h"
#include "undergen_smooth_node.h"
#include "undergen_vox_stamp_node.h"
#include "undergen_surface_sampler_node.h"
#include "undergen_point_filter_node.h"
#include "undergen_mesher_node.h"
#include "undergen_scene_spawner_node.h"
#include "undergen_grammar_node.h"
#include "undergen_material_stamper_node.h"
#include "undergen_liquid_flood_node.h"
#include "undergen_grid_warp_node.h"
#include "zone_material_entry.h"
#include "vox_spawn_entry.h"
#include "vox_material_entry.h"

#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/core/defs.hpp>
#include <gdextension_interface.h>
#include <godot_cpp/godot.hpp>

using namespace godot;

void initialize_density_grid_module(ModuleInitializationLevel p_level) {
    if (p_level != MODULE_INITIALIZATION_LEVEL_SCENE) {
        return;
    }

    // --- Core classes ---
    ClassDB::register_class<DensityGrid>();
    ClassDB::register_class<MCChunk>();

    // --- New Pipeline Infrastructure ---
    ClassDB::register_class<UnderGenPointSet>();
    ClassDB::register_class<UnderGenNode>();
    ClassDB::register_class<UnderGenPipeline>();
    ClassDB::register_class<UnderGenWorld3D>();

    // --- Grammar & Pipeline Spec (code/JSON builder API) ---
    ClassDB::register_class<LevelGrammarRoomTypeSpec>();
    ClassDB::register_class<LevelGrammarRuleSpec>();
    ClassDB::register_class<LevelGrammarSpec>();

    // --- Concrete Pipeline Nodes ---
    ClassDB::register_class<UnderGenBSPPlacerNode>();
    ClassDB::register_class<UnderGenAStarCarverNode>();
    ClassDB::register_class<UnderGenBezierCarverNode>();
    ClassDB::register_class<UnderGenNoiseNode>();
    ClassDB::register_class<UnderGenSmoothNode>();
    ClassDB::register_class<UnderGenVoxStampNode>();
    ClassDB::register_class<UnderGenSurfaceSamplerNode>();
    ClassDB::register_class<UnderGenPointFilterNode>();
    ClassDB::register_class<UnderGenMesherNode>();
    ClassDB::register_class<UnderGenSceneSpawnerNode>();
    ClassDB::register_class<UnderGenGrammarNode>();
    ClassDB::register_class<UnderGenMaterialStamperNode>();
    ClassDB::register_class<UnderGenLiquidFloodNode>();
    ClassDB::register_class<UnderGenGridWarpNode>();
    ClassDB::register_class<ZoneMaterialEntry>();
    ClassDB::register_class<VoxSpawnEntry>();
    ClassDB::register_class<VoxMaterialEntry>();
}

void uninitialize_density_grid_module(ModuleInitializationLevel p_level) {
    if (p_level != MODULE_INITIALIZATION_LEVEL_SCENE) {
        return;
    }
}

extern "C" {
GDExtensionBool GDE_EXPORT gdextension_entry_point(
    GDExtensionInterfaceGetProcAddress p_get_proc_address,
    GDExtensionClassLibraryPtr p_library,
    GDExtensionInitialization *r_initialization)
{
    godot::GDExtensionBinding::InitObject init_obj(p_get_proc_address, p_library, r_initialization);
    init_obj.register_initializer(initialize_density_grid_module);
    init_obj.register_terminator(uninitialize_density_grid_module);
    init_obj.set_minimum_library_initialization_level(MODULE_INITIALIZATION_LEVEL_SCENE);
    return init_obj.init();
}
} // extern "C"