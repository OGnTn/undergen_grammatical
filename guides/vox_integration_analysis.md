# Analysis: Integrating MagicaVoxel (.vox) Structures

## Overview
This document analyzes how to integrate premade MagicaVoxel (`.vox`) structures into the existing `Logical Grammar -> DensityGrid -> Spawning` pipeline.

The goal is to allow Game Designers to inspect and modify rooms visually in MagicaVoxel, defining both the **Geometry** (walls, floors) and **Gameplay Elements** (spawn points, chests) using voxel colors, and have the procedural engine respect these definitions.

## Workflow (Game Designer Perspective)
1.  **Creation**: You open MagicaVoxel and model a `CryptRoom`.
    *   **Geometry**: You paint the walls and floor.
    *   **Spawning**: You use specific palette colors to mark locations.
        *   *Example*: Green Pixel = "Zombie Spawn". Gold Pixel = "Chest".
2.  **Definition**: In Godot's Grammar Editor, you create a symbol `Crypt`.
    *   You assign `vox_file = "res://assets/vox/crypt_room.vox"`.
    *   You map Palette Indices to Gameplay concepts (e.g., Index 255 -> Zombie).
3.  **Generation**: You run the generator. The system logically places a `Crypt` node, stamps the voxel geometry into the terrain, and instantiates Zombies at the exact pixel locations you painted.

## Technical Integration Points

### 1. Logical Grammar (The "Brain")
**File**: `GrammarGenerator.gd` / `LevelGrammarResource`
*   **Current State**: Nodes have `min_size`, `max_size`, `type`.
*   **Change**: Add a property `vox_path` (String) to the Grammar Symbol definitions.
*   **Flow**: When `GrammarGenerator` creates a Node for `Crypt`, it passes the `vox_path` string in the node dictionary returned to `World.gd`.

### 2. Density Grid (The "Canvas")
**File**: `LevelDensityGrid.cpp` (C++)
*   **Current State**: `generate_from_graph` calculates a bounding box and runs `RoomGenerator` to carve generic shapes (boxes).
*   **Change**: 
    *   **Include Library**: Use the existing `src/ogt_vox.h` library (Single-header C++ library).
    *   **Update `ResolvedRoom`**: Add `vox_path` member.
    *   **Vox Cache**: Implement a static map `std::map<String, const ogt_vox_scene*>` to cache loaded scenes and avoid disk I/O spam.
    *   **Loading**: Use `ogt_vox_read_scene(buffer, size)` to parse `.vox` files.
    *   **Geometry Stamping**: 
        *   **Condition**: Check if `vox_path` is not empty and file exists.
        *   **If Vox Exists**:
            *   Iterate `scene->instances`.
            *   For each instance, access its `model` and `voxel_data`.
            *   **Palette Check**: 
                *   Get `palette_index` of voxel.
                *   **Case A: Spawn Point**: If index is in `vox_spawn_map`, record spawn and SKIP density write (treat as air).
                *   **Case B: Solid Block**: 
                    *   Set Density to `SOLID` (0.0).
                    *   **Material ID**: Lookup `vox_material_map` (Index -> Material ID). Apply this Material ID to the grid cell (requires `set_material_id` support in DensityGrid).
            *   **Empty Voxel** (0): Set Density to `OPEN` (1.0).
        *   **Fallback (If No Vox)**:
            *   Call the existing carving logic.

### 3. Spawning (The "Inhabitants")
**File**: `populator.gd`
*   **Current State**: Uses Poisson Disk Sampling on surface normals.
*   **Change**: Implement **Hybrid Spawning**.
    1.  **Fixed Spawns (Vox)**:
        *   Call `LevelDensityGrid.get_vox_spawns()`. 
        *   Instantiate these objects first.
    2.  **Procedural Spawns (Surface)**:
        *   Run the existing Poisson Disk Sampler.
        *   *Conflict Avoidance*: Ensure Poisson candidates don't overlap with Fixed Spawns (optional check, or just rely on density checks).
    *   **Palette Allocation**:
        *   Reserve specific indices (e.g., 240-255) for Spawns.
        *   Use indices 1-239 for Materials (Walls, Floors, Trim).
        *   Define this range in `LevelResource` to avoid overlap.

## Implementation Roadmap

### Phase 1: Data Pipeline
1.  **Grammar**: Add `vox_file` export to `GrammarSymbol` resource.
2.  **Transfer**: Ensure `GrammarGenerator` passes `vox_file` to `ResolvedRoom` in C++.

### Phase 2: Geometry Stamping (using `ogt_vox`)
1.  **Integration**: 
    *   Include `ogt_vox.h` in `level_density_grid.cpp`.
    *   Create a `load_vox(path)` helper that reads the file into a buffer and calls `ogt_vox_read_scene`.
    *   Manage memory: `ogt_vox_destroy_scene` when the grid is destroyed or cache cleared.
2.  **RoomGenerator**: Add `stamp_vox(LevelDensityGrid* grid, ResolvedRoom room, const ogt_vox_scene* scene)` method.
    *   Map Voxel Space -> Grid Space, applying Room Position offset.
    *   Handle Scale: Assume 1 Voxel = 1 Godot Unit for now (matches `cube_size`).

### Phase 3: Metadata & Spawning
1.  **Palette Config**: Add `vox_palette_map` to `LevelResource` (Dictionary: Index -> SpawnableObject).
2.  **Extraction**: During `stamp_vox`, if index matches a Key in `vox_palette_map`:
    *   Don't write density (treat as air).
    *   Add to `spawn_points` list.
3.  **Populator**: Read `spawn_points` and instantiate scenes.

## Considerations
*   **Scaling**: `.vox` coordinates are integers. Ensure `1 voxel` = `1 grid cell` (usually 1 meter) or define a scale factor.
*   **Rotation**: If Grammar rotates the room (e.g. 90 degrees), the Vox data must be rotated before stamping.
*   **Performance**: Reading `.vox` files from disk is fast, but avoid re-reading the same file 100 times. Cache the parsed voxel data in a static map in C++.
