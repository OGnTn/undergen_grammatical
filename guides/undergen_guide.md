# UnderGen System Architecture Guide

## 1. High-Level Workflow
The generation pipeline follows a strict **Logic -> Physics -> Gameplay** flow.

```mermaid
graph TD
    A[Grammar Resource] -->|Rules & Logic| B[GrammarGenerator (GDScript)]
    B -->|Graph: Nodes & Edges| C[LevelDensityGrid (C++)]
    C -->|Voxel Data| D[Mesh Generation]
    C -->|Zone IDs| E[Gameplay Spawner]
```

## 2. The Logic Layer: Grammar
This is where you **architect** the level's flow. You do not think about coordinates here, only relationships.
*   **Symbols**: Abstract tags (e.g., `Start`, `Boss`, `Treasure`).
*   **Rules**: Recipes for replacing symbols.
    *   *Example*: `Hallway` -> `HallwaySegment` + `Enemy` + `HallwaySegment`.
*   **State**: Variables (e.g., `keys`, `danger_level`) that you can:
    *   **Check (Condition)**: `keys >= 3`
    *   **Modify (Action)**: `{"keys": 1}`

## 3. The Physical Layer: Density Grid
Once the logical graph is built, it is sent to the C++ engine. This layer creates the **Volume**.
*   **Room Generator**: Takes nodes and places them as rooms (finding non-overlapping spots).
    *   *Control*: `min_room_size`, `max_room_size` in the resource.
*   **Path Carver**: Connects the rooms based on the graph edges.
    *   *Control*: `path_radius`, `dungeon_mode` (straight vs bezier).
*   **Zone System**:
    *   Each node has a **Zone Name** (from its Symbol).
    *   The grid stores this name in every voxel belonging to that room.
    *   *Result*: You know exactly which room is where.

## 4. The Gameplay Layer: Spawning
This is where you turn geometry into a game. You use a weighted **Poisson Disk sampler** (`Populator`) that places objects on valid surfaces.

### How to Place Objects
Instead of manual placement coordinates, you define rules in your `WorldSpawnableObject` resources:

1.  **Define the Asset**: Create a resource for your "Chest".
2.  **Set Zone Rules**:
    *   Find the `Required Zones` array in the inspector.
    *   Add `"ChestRoomNode"`.
3.  **The System Does the Rest**:
    *   The `Populator` scans the level surface.
    *   It checks `density_grid.get_zone_at(position)`.
    *   If the zone matches your requirement, the chest has a chance to spawn.

*This allows for organic placement (e.g., "Ferns only in 'Jungle' zone") without hardcoded positions.*

## 5. What You Can Architect
| Feature | Controlled By | Description |
| :--- | :--- | :--- |
| **Pacing / Flow** | Grammar Rules | Ensure `Key` appears before `Door`. Control branching depth. |
| **Room Shapes** | LevelDensityGrid | `min/max_room_size` controls variance. `dungeon_mode` toggles caves vs halls. |
| **Corridor Style** | Edge Types | Tag edges as `bridge` or `tunnel` in editor to change connectivity feel. |
| **Difficulty** | Logic Variables | Increment `difficulty` var each room, spawn harder enemies later. |

## 6. Full Example Pipeline
1.  **Editor**: Create `LevelGrammar` resource.
    *   Rule: `Start` -> `Entry` -- `Hallway`.
    *   Rule: `Hallway` -> `Room` -- `Hallway` (Recursive).
    *   Rule: `Room` -> `Treasure` (Low Prob).
    *   Rule: `Room` -> `MobArena` (High Prob).
2.  **Game**:
    *   `world.gd` calls `generator.generate()`.
    *   System places `Entry`, `MobArena`, `MobArena`, `Treasure` physically.
    *   System carves paths between them.
    *   `Populator` detects "MobArena" zone and spawns enemies.
    *   `Populator` detects "Treasure" zone and spawns loot.
