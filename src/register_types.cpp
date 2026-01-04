// src/register_types.cpp
#include "register_types.h"

#include "density_grid.h"
#include "level_density_grid.h"
#include "mc_chunk.h"
#include "chunk_manager.h"

#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/core/defs.hpp>
#include <gdextension_interface.h> // Correct header for GDExtensionInterface
#include <godot_cpp/godot.hpp>

using namespace godot;

// Module initialization callback function
void initialize_density_grid_module(ModuleInitializationLevel p_level) {
    if (p_level != MODULE_INITIALIZATION_LEVEL_SCENE) {
        return; // Only initialize when the scene level is reached
    }

    // Register custom classes
    ClassDB::register_class<DensityGrid>();
    ClassDB::register_class<LevelDensityGrid>();
    ClassDB::register_class<MCChunk>();
    ClassDB::register_class<ChunkManager>();
    // You could register singletons or editor plugins here if needed
}

// Module termination callback function
void uninitialize_density_grid_module(ModuleInitializationLevel p_level) {
    if (p_level != MODULE_INITIALIZATION_LEVEL_SCENE) {
        return;
    }
    // Uninitialization logic if necessary (usually not needed for simple classes)
}

// GDExtension entry point definition
// Tells Godot where to find the initialization and termination functions
extern "C" {
// Use GDExtensionBool to match the expected return type
GDExtensionBool GDE_EXPORT gdextension_entry_point(
    GDExtensionInterfaceGetProcAddress p_get_proc_address,
    GDExtensionClassLibraryPtr p_library,
    GDExtensionInitialization *r_initialization)
{
    godot::GDExtensionBinding::InitObject init_obj(p_get_proc_address, p_library, r_initialization);

    // Set up the initialization and termination function pointers
    init_obj.register_initializer(initialize_density_grid_module);
    init_obj.register_terminator(uninitialize_density_grid_module);
    init_obj.set_minimum_library_initialization_level(MODULE_INITIALIZATION_LEVEL_SCENE);

    return init_obj.init(); // Finalize the initialization object
}
} // extern "C"