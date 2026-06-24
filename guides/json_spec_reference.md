# UnderGen JSON Specification v1.0

Complete reference for defining **grammars** and **pipelines** as JSON.

---

## Table of Contents

1. [Grammar Spec (`*.grammar.json`)](#1-grammar-spec-grammarjson)
   - [Top-Level Object](#11-top-level-object)
   - [Room Type Object](#12-room-type-object)
   - [Rule Object](#13-rule-object)
   - [RHS Node Object](#14-rhs-node-object)
   - [RHS Edge Object](#15-rhs-edge-object)
   - [Action Object](#16-action-object)
   - [Condition Operators](#17-condition-operators)
2. [Pipeline Spec (`*.pipeline.json`)](#2-pipeline-spec-pipelinejson)
   - [Top-Level Object](#21-top-level-object)
   - [Pipeline Node Object](#22-pipeline-node-object)
   - [Connection Object](#23-connection-object)
3. [Node Type Reference](#3-node-type-reference)
   - [Port Type Legend](#31-port-type-legend)
   - [UnderGenGrammarNode — Grammar Expander](#32-undergengrammarnode)
   - [UnderGenBSPPlacerNode — BSP Room Placer](#33-undergenbspplacernode)
   - [UnderGenAStarCarverNode — A* Corridor Carver](#34-undergenastarcarvernode)
   - [UnderGenBezierCarverNode — Bezier Tunnel Carver](#35-undergenbeziercarvernode)
   - [UnderGenNoiseNode — 3D Noise Blend](#36-undergennoisenode)
   - [UnderGenSmoothNode — Smooth Filter](#37-undergensmoothnode)
   - [UnderGenVoxStampNode — Vox Stamp](#38-undergenvoxstampnode)
   - [UnderGenMaterialStamperNode — Material Stamper](#39-undergenmaterialstampernode)
   - [UnderGenSurfaceSamplerNode — Surface Sampler](#310-undergensurfacesamplernode)
   - [UnderGenPointFilterNode — Point Filter](#311-undergenpointfilternode)
   - [UnderGenMesherNode — Mesher (Terminal)](#312-undergenmeshernode)
   - [UnderGenSceneSpawnerNode — Scene Spawner (Terminal)](#313-undergenscenespawnernode)
4. [Full Examples](#4-full-examples)
5. [API Usage](#5-api-usage)

---

## 1. Grammar Spec (`*.grammar.json`)

Defines a **level grammar**: an axiom, a palette of room symbols, and expansion rules.

### 1.1 Top-Level Object

```jsonc
{
  // ── Required ──────────────────────────────────────────────
  "axiom": "Start",                     // String — the initial symbol

  // ── Optional ──────────────────────────────────────────────
  "state_variables": ["keys", "level"], // Array<String> — global state vars

  // ── Required ──────────────────────────────────────────────
  "room_types": [ /* ... */ ],          // Array<RoomType> — palette of symbols
  "rules":       [ /* ... */ ]          // Array<Rule>     — expansion rules
}
```

| Key | Type | Default | Description |
|:----|:-----|:--------|:------------|
| `axiom` | `string` | — | The starting symbol. One root node with this symbol is created before iteration 1. |
| `state_variables` | `string[]` | `[]` | Names of global state variables tracked during expansion (e.g. `"keys"`, `"difficulty"`). Used in rule conditions and actions. |
| `room_types` | `RoomType[]` | — | The palette. Every symbol referenced in rules must appear here. Determines spatial constraints and visual identity. |
| `rules` | `Rule[]` | — | Expansion rules. Applied iteratively: at each iteration, every node that matches a rule's `lhs_symbol` may be replaced by that rule's RHS subgraph. |

---

### 1.2 Room Type Object

```jsonc
{
  "symbol":   "Entry",                  // String — unique symbol identifier
  "color":    [0.0, 1.0, 1.0, 1.0],    // Color as [R, G, B, A] floats 0..1
  "weight":   1.0,                      // float — weight for placement solver

  // Spatial bounds (integer voxels, Vector3i)
  "min_size": [4, 3, 4],               // [x, y, z] — smallest allowed room
  "max_size": [8, 5, 8],               // [x, y, z] — largest allowed room

  "vox_path": "res://vox/entry.vox"    // String — optional .vox file to stamp
}
```

| Key | Type | Default | Description |
|:----|:-----|:--------|:------------|
| `symbol` | `string` | — | Unique name. Rules reference it via `lhs_symbol`. Terminal nodes inherit their symbol as the zone name. |
| `color` | `[float,float,float,float]` | `[0.4,0.6,0.9,1]` | Editor display colour. RGBA, each 0..1. |
| `weight` | `float` | `1.0` | Weighting bias for the BSP room placer. |
| `min_size` | `[int,int,int]` | `[5,3,5]` | Minimum room dimensions in voxels `[x, y, z]`. |
| `max_size` | `[int,int,int]` | `[10,6,10]` | Maximum room dimensions in voxels `[x, y, z]`. Actual size is randomly chosen between these bounds per room. |
| `vox_path` | `string` | `""` | Path to a `.vox` (MagicaVoxel) file stamped into the room. Empty = no vox. |

---

### 1.3 Rule Object

```jsonc
{
  "rule_name":    "Main Split",         // String — human-readable label
  "lhs_symbol":   "Start",              // String — symbol this rule replaces
  "probability":  1.0,                  // float — weight for random selection

  "entry_node_id": "entry",             // String — which RHS node links to parent's incoming edges
  "exit_node_id":  "boss",              // String — which RHS node links to parent's outgoing edges

  // ── Condition (all optional — omit entire block for unconditional) ──
  "condition_var": "keys",              // String — state variable to check
  "condition_op":  ">=",                // String — operator (see §1.7)
  "condition_val": 1,                   // float  — threshold value

  // ── Actions ───────────────────────────────────────────────
  "actions": [                          // Action[] — state mutations on rule fire
    { "var": "keys", "delta": 1 }
  ],

  // ── RHS subgraph ──────────────────────────────────────────
  "rhs_nodes": [ /* ... */ ],           // RHSNode[] — nodes in replacement subgraph
  "rhs_edges": [ /* ... */ ]            // RHSEdge[] — edges in replacement subgraph
}
```

| Key | Type | Default | Description |
|:----|:-----|:--------|:------------|
| `rule_name` | `string` | `"New Rule"` | Display name. |
| `lhs_symbol` | `string` | — | The symbol to replace. Must match one in `room_types`. |
| `probability` | `float` | `1.0` | Weight when multiple rules match the same LHS. Higher = more likely. |
| `entry_node_id` | `string` | `""` | ID of the RHS node that receives edges previously connected to the replaced node. If empty, first RHS node is used. |
| `exit_node_id` | `string` | `""` | ID of the RHS node from which outgoing edges leave. If empty, same as entry. |
| `condition_var` | `string` | `""` | State variable to test. Empty = always fires. |
| `condition_op` | `string` | `"<"` | Comparison operator. See §1.7. |
| `condition_val` | `float` | `0` | Comparison threshold. |
| `actions` | `Action[]` | `[]` | State mutations applied when this rule fires. See §1.6. |
| `rhs_nodes` | `RHSNode[]` | — | Nodes of the replacement subgraph. See §1.4. |
| `rhs_edges` | `RHSEdge[]` | — | Edges of the replacement subgraph. See §1.5. |

> **How rules work**: At each iteration, for each node in the graph, all rules whose `lhs_symbol` matches are checked. Among those whose condition passes (or is empty), one is chosen randomly weighted by `probability`. The node is removed and replaced by the rule's RHS subgraph. External edges are rewired: incoming edges → `entry_node_id`, outgoing edges ← `exit_node_id`. If no rule matches a node's symbol, it becomes a **terminal leaf** — it stays in the graph and eventually becomes a placed room.

---

### 1.4 RHS Node Object

```jsonc
{
  "id":     "entry",                    // String — local ID (scoped within this rule)
  "symbol": "Entry",                    // String — must appear in room_types palette

  "constraints": {                      // Dict — optional spatial constraints
    "fixed_pos":  [0, 0, 0],           // [x, y, z] float — lock to a grid-relative position
    "relative_to": "hall"              // String — constrain position relative to another node in the SAME rule
  }
}
```

| Key | Type | Default | Description |
|:----|:-----|:--------|:------------|
| `id` | `string` | — | Unique within this rule. Referenced by `rhs_edges`, `entry_node_id`, `exit_node_id`, and `constraints.relative_to`. |
| `symbol` | `string` | — | Room symbol. Must exist in `room_types`. If this symbol has no matching rules, it becomes terminal. |
| `constraints` | `object` | `{}` | Placement constraints for the BSP placer and path carver. |
| `constraints.fixed_pos` | `[float,float,float]` | — | Lock the room at this position relative to the grid center. |
| `constraints.relative_to` | `string` | — | Place this node relative to another node (by local ID) in the same RHS. |

---

### 1.5 RHS Edge Object

```jsonc
{
  "from": "entry",                      // String — local ID of source RHS node
  "to":   "hall",                       // String — local ID of target RHS node
  "type": "corridor"                    // String — edge type for path carving
}
```

| Key | Type | Default | Description |
|:----|:-----|:--------|:------------|
| `from` | `string` | — | Source RHS node `id`. |
| `to` | `string` | — | Target RHS node `id`. |
| `type` | `string` | `"corridor"` | Edge class. Can be `"corridor"`, `"bridge"`, `"tunnel"`, or custom. Affects carving style and zone. |

---

### 1.6 Action Object

```jsonc
{ "var": "keys", "delta": 1 }          // Adds 1 to the "keys" state variable
```

| Key | Type | Description |
|:----|:-----|:------------|
| `var` | `string` | State variable name. |
| `delta` | `float` | Value to add (can be negative). Use `0` in conjunction with a non-numeric value to set it directly. |

---

### 1.7 Condition Operators

| Operator | Meaning |
|:---------|:--------|
| `"<"` | less than |
| `">"` | greater than |
| `"<="` | less than or equal |
| `">="` | greater than or equal |
| `"=="` | equal (within 1e-6 tolerance) |
| `"!="` | not equal |

If `condition_var` is empty (`""`), the rule always fires (no condition check).

---

## 2. Pipeline Spec (`*.pipeline.json`)

Defines a **processing pipeline**: a directed acyclic graph (DAG) of generation nodes.

### 2.1 Top-Level Object

```jsonc
{
  "nodes":       [ /* PipelineNode[] */ ],    // generation steps
  "connections": [ /* Connection[]   */ ]     // wires between steps
}
```

| Key | Type | Description |
|:----|:-----|:------------|
| `nodes` | `PipelineNode[]` | Ordered list of processing nodes. Each has a `type`, `name`, and `properties`. See §2.2. |
| `connections` | `Connection[]` | Directed edges connecting node output ports to input ports. See §2.3. |

---

### 2.2 Pipeline Node Object

```jsonc
{
  "name": "Grammar",                    // String — unique name within pipeline
  "type": "UnderGenGrammarNode",       // String — C++ class name (see §3 for all types)
  "pos":  [100, 200],                   // [float, float] — editor canvas position (optional, for visual layout)

  "properties": {                       // Dict — node-specific settings (see §3 per type)
    "grammar_resource_path": "res://grammars/dungeon.tres",
    "iterations": 5,
    "max_nodes": 150
  }
}
```

| Key | Type | Default | Description |
|:----|:-----|:--------|:------------|
| `name` | `string` | — | Unique identifier. Connections reference this name. |
| `type` | `string` | — | C++ class name. Must be one of the types listed in §3. |
| `pos` | `[float,float]` | `[0,0]` | Position on the editor canvas. Purely cosmetic. |
| `properties` | `object` | `{}` | Key-value map of the node's exported properties. See §3 for the schema of each type. |

---

### 2.3 Connection Object

```jsonc
{
  "from_node": "Grammar",              // String — source node name
  "from_port": 0,                       // int    — source output port index
  "to_node":   "BSP_Placer",           // String — target node name
  "to_port":    1                       // int    — target input port index
}
```

| Key | Type | Description |
|:----|:-----|:------------|
| `from_node` | `string` | Node name of the source. |
| `from_port` | `int` | Output port index on the source node (see §3 per type). |
| `to_node` | `string` | Node name of the destination. |
| `to_port` | `int` | Input port index on the destination node (see §3 per type). |

---

## 3. Node Type Reference

Every pipeline node has a set of **typed ports** and **configurable properties**.

### 3.1 Port Type Legend

| Type ID | Name | Description | Voxel color |
|:--------|:-----|:------------|:------------|
| `0` | **Seed** | Integer seed for RNG. Used as initial pipeline input. | White |
| `1` | **Logical Graph** | `Dictionary { "nodes": Array, "edges": Array }` — the grammar expansion output. | Blue |
| `2` | **Gen Context** | `Dictionary { "grid": DensityGrid, "rooms": Array, "edges": Array, ... }` — the shared generation state flowing through the pipeline. | Green |
| `3` | **Point Set** | `UnderGenPointSet` — collection of surface sample points with attributes. | Yellow |

---

### 3.2 UnderGenGrammarNode

**Editor Label**: Grammar Expander
**Role**: Runs iterative grammar expansion to produce a logical graph.

| Port | Direction | Type | Description |
|:-----|:----------|:-----|:------------|
| 0 | **in** | `0` Seed | Integer — seeds the RNG. |
| 0 | **out** | `1` Logical Graph | The expansion result. |

| Property | Type | Default | Range | Description |
|:---------|:-----|:--------|:------|:------------|
| `grammar_resource_path` | `string` | `""` | `*.tres,*.res` | Path to the `LevelGrammarResource` or `LevelGrammarSpec` `.tres` file. |
| `iterations` | `int` | `4` | `1..32` | Number of expansion passes. Higher = deeper nesting, more rooms. |
| `max_nodes` | `int` | `100` | `1..1000` | Hard cap on total nodes. Expansion stops when reached. |

---

### 3.3 UnderGenBSPPlacerNode

**Editor Label**: BSP Room Placer
**Role**: Places rooms from the logical graph into the voxel density grid using Binary Space Partitioning.

| Port | Direction | Type | Description |
|:-----|:----------|:-----|:------------|
| 0 | **in** | `0` Seed | Integer — seeds the RNG. |
| 1 | **in** | `1` Logical Graph | The grammar output. |
| 0 | **out** | `2` Gen Context | Grid + placed rooms + edges. |

| Property | Type | Default | Range | Description |
|:---------|:-----|:--------|:------|:------------|
| `grid_size_x` | `int` | `48` | — | Grid width in voxels. |
| `grid_size_y` | `int` | `16` | — | Grid height in voxels. |
| `grid_size_z` | `int` | `48` | — | Grid depth in voxels. |

---

### 3.4 UnderGenAStarCarverNode

**Editor Label**: A* Carver
**Role**: Carves straight/direct corridors between rooms using 3D A\* pathfinding.

| Port | Direction | Type | Description |
|:-----|:----------|:-----|:------------|
| 0 | **in** | `2` Gen Context | Grid with placed rooms. |
| 0 | **out** | `2` Gen Context | Grid with carved paths. |

| Property | Type | Default | Range | Description |
|:---------|:-----|:--------|:------|:------------|
| `path_brush_min_radius` | `int` | `2` | — | Minimum corridor radius in voxels. |
| `path_brush_max_radius` | `int` | `4` | — | Maximum corridor radius in voxels. Random between min/max per segment. |
| `use_square_brush` | `bool` | `false` | — | If `true`, carves square-profile corridors instead of round. |
| `vertical_movement_cost_multiplier` | `float` | `2.0` | — | A\* cost multiplier for vertical movement. Higher = paths stay more horizontal. |

---

### 3.5 UnderGenBezierCarverNode

**Editor Label**: Bezier Carver
**Role**: Carves organic, winding tunnels using bezier curves with configurable wobble and cave shaping.

| Port | Direction | Type | Description |
|:-----|:----------|:-----|:------------|
| 0 | **in** | `2` Gen Context | Grid with placed rooms. |
| 0 | **out** | `2` Gen Context | Grid with carved tunnels. |

#### Brush Group (prefix: `path_brush_`)

| Property | Type | Default | Range | Description |
|:---------|:-----|:--------|:------|:------------|
| `path_brush_min_radius` | `int` | `2` | — | Minimum corridor radius. |
| `path_brush_max_radius` | `int` | `4` | — | Maximum corridor radius. |
| `path_brush_square` | `bool` | `false` | — | Square brush profile. |

#### Bezier Group (prefix: `path_`)

| Property | Type | Default | Range | Description |
|:---------|:-----|:--------|:------|:------------|
| `path_segments` | `int` | `1` | `1..10` | Number of bezier segments per path. More = more control points and snaking. |
| `path_bend_factor` | `float` | `0.4` | `0..2` | How far control points can deviate from the straight line. Higher = wilder curves. |
| `path_wobble_magnitude` | `float` | `0` | `0..10` | Amplitude of noise-based wobble along the path. |
| `path_wobble_frequency` | `float` | `0.2` | `0..1` | Frequency of the wobble noise. Higher = more tight wiggles. |
| `path_connect_from_ground_level` | `bool` | `false` | — | If `true`, connections start/end at room floor level instead of room center. |

#### Cave Shape Group (prefix: `cave_`)

| Property | Type | Default | Range | Description |
|:---------|:-----|:--------|:------|:------------|
| `cave_ruggedness` | `float` | `1.0` | `0..5` | Overall roughness of cave walls. 0 = smooth tube. |
| `cave_floor_ruggedness` | `float` | `0` | `0..5` | Extra roughness applied to the floor surface. |
| `cave_ceiling_ruggedness` | `float` | `0` | `0..5` | Extra roughness applied to the ceiling surface. |
| `cave_width_noise` | `float` | `0` | `0..5` | Variation in corridor width along the path. |
| `floor_flattening` | `float` | `0` | `0..1` | Bias toward flat floors. 1.0 = completely flat. |
| `overhang_openness` | `float` | `0` | `0..2` | Tendency for overhanging rock formations. 0 = none. |

---

### 3.6 UnderGenNoiseNode

**Editor Label**: 3D Noise Blend
**Role**: Blends 3D Perlin/Simplex noise into the density grid for organic cave feel.

| Port | Direction | Type | Description |
|:-----|:----------|:-----|:------------|
| 0 | **in** | `2` Gen Context | Grid to noise-blend. |
| 0 | **out** | `2` Gen Context | Grid with noise applied. |

| Property | Type | Default | Range | Description |
|:---------|:-----|:--------|:------|:------------|
| `noise_scale` | `float` | `50.0` | — | World-space scale of the noise. Larger = broader features. |
| `noise_intensity` | `float` | `0.5` | — | Blend strength. 0 = no effect, 1 = full amplitude. |
| `noise_frequency` | `float` | `0.02` | — | Base frequency of the noise generator. |
| `noise_seed` | `int` | `1337` | — | Noise RNG seed. |
| `noise_generator` | `FastNoiseLite` | *(auto)* | — | Godot `FastNoiseLite` resource. If omitted, auto-created. Set path as `"res://noise/my_noise.tres"` in JSON. |

---

### 3.7 UnderGenSmoothNode

**Editor Label**: Smooth Filter
**Role**: Applies a separable 3D box-blur low-pass filter to the density grid.

| Port | Direction | Type | Description |
|:-----|:----------|:-----|:------------|
| 0 | **in** | `2` Gen Context | Grid to smooth. |
| 0 | **out** | `2` Gen Context | Smoothed grid. |

| Property | Type | Default | Range | Description |
|:---------|:-----|:--------|:------|:------------|
| `smoothing_strength` | `int` | `1` | — | Box-blur radius in voxels. 1 = 3×3×3 kernel, 2 = 5×5×5. Higher = softer geometry. |

---

### 3.8 UnderGenVoxStampNode

**Editor Label**: Vox Stamper
**Role**: Stamps `.vox` (MagicaVoxel) models into rooms and extracts spawn-point markers.

| Port | Direction | Type | Description |
|:-----|:----------|:-----|:------------|
| 0 | **in** | `2` Gen Context | Grid with placed rooms. |
| 0 | **out** | `2` Gen Context | Grid with voxel stamps applied. |

| Property | Type | Default | Description |
|:---------|:-----|:--------|:------------|
| `vox_spawn_map` | `Dict<int → string>` | `{}` | Map palette index → spawn marker type. e.g. `{ "3": "enemy_spawn", "5": "loot_spawn" }`. |
| `vox_material_map` | `Dict<int → int>` | `{}` | Map palette index → material ID. e.g. `{ "1": 0, "2": 1 }`. |
| `vox_inverse_density` | `bool` | `false` | If `true`, voxels carve (set to 0) instead of filling (set to 1). |

---

### 3.9 UnderGenMaterialStamperNode

**Editor Label**: Material Stamper
**Role**: Assigns material IDs to voxels based on zone name → material ID mapping.

| Port | Direction | Type | Description |
|:-----|:----------|:-----|:------------|
| 0 | **in** | `2` Gen Context | Grid to stamp materials on. |
| 0 | **out** | `2` Gen Context | Grid with material IDs set. |

| Property | Type | Default | Range | Description |
|:---------|:-----|:--------|:------|:------------|
| `zone_material_map` | `ZoneMaterialEntry[]` | `[]` | — | Array of `{ "zone_name": string, "material_id": int }` entries. See below. |
| `default_material_id` | `int` | `0` | `0..255` | Material ID for any voxel whose zone is not in the map. |

**ZoneMaterialEntry** (each element of `zone_material_map`):
```jsonc
{ "zone_name": "bossroom", "material_id": 3 }
```

---

### 3.10 UnderGenSurfaceSamplerNode

**Editor Label**: Surface Sampler
**Role**: Samples voxel surfaces (floors, ceilings, walls) into a `PointSet` for object spawning.

| Port | Direction | Type | Description |
|:-----|:----------|:-----|:------------|
| 0 | **in** | `2` Gen Context | Grid to sample. |
| 0 | **out** | `2` Gen Context | Pass-through. |
| 1 | **out** | `3` Point Set | Sampled surface points. |

| Property | Type | Default | Range | Description |
|:---------|:-----|:--------|:------|:------------|
| `surface_type` | `int` (enum) | `0` | `0=Floor, 1=Ceiling, 2=Wall, 3=All` | Which surface normals to sample. |
| `slope_threshold` | `float` | `0.6` | — | Dot-product threshold for floor vs wall classification. |
| `voxel_size` | `float` | `1.0` | — | World size of one voxel (for world-space position calculation). |
| `zone_filter` | `string` | `""` | Comma-separated | Only sample voxels in these zones. Empty = all zones. e.g. `"bossroom,corridor"`. |
| `zone_match_mode` | `int` (enum) | `0` | `0=Exact, 1=Prefix` | `Exact`: full zone name match. `Prefix`: zone name starts with filter. |

---

### 3.11 UnderGenPointFilterNode

**Editor Label**: Point Filter
**Role**: Filters a PointSet by zone, slope, density threshold, and minimum spacing. Used upstream of `UnderGenSceneSpawnerNode` to route points by zone into separate spawners.

| Port | Direction | Type | Description |
|:-----|:----------|:-----|:------------|
| 0 | **in** | `3` Point Set | Points to filter. |
| 0 | **out** | `3` Point Set | Filtered points. |

| Property | Type | Default | Range | Description |
|:---------|:-----|:--------|:------|:------------|
| `required_zone_name` | `string` | `""` | Comma-separated | Keep or exclude points in these zones (see `zone_match_mode`). Empty = all. |
| `zone_match_mode` | `int` (enum) | `0` | `0=Exact, 1=Prefix, 2=Exclude` | How to match zone names. `Exact`: full name match. `Prefix`: zone name starts with filter. `Exclude`: keep points whose zone matches **none** of the names. |
| `min_slope` | `float` | `0` | `0..1` | Minimum slope (0=flat, 1=vertical). |
| `max_slope` | `float` | `1` | `0..1` | Maximum slope. |
| `min_density` | `float` | `0` | — | Points with `density` below this are discarded. |
| `min_spacing` | `float` | `0` | World units | Minimum distance between kept points. 0 = no spacing check. |

---

### 3.12 UnderGenMesherNode

**Editor Label**: Marching Cubes Mesher
**Role**: **Terminal node.** Spawns `MCChunk` mesh nodes from the density grid. Requires `UnderGenWorld3D` as parent.

| Port | Direction | Type | Description |
|:-----|:----------|:-----|:------------|
| 0 | **in** | `2` Gen Context | Final grid to mesh. |
| *(no output)* | | | End of geometry pipeline. |

| Property | Type | Default | Range | Description |
|:---------|:-----|:--------|:------|:------------|
| `chunk_size` | `int` | `16` | — | Voxels per chunk (cubed). 16 = 16×16×16 chunks. |
| `voxel_size` | `float` | `1.0` | — | World size of one voxel in units. |
| `generate_collision` | `bool` | `true` | — | Create collision shapes for the generated meshes. |
| `generate_occluder` | `bool` | `false` | — | Create occlusion culling shapes. |
| `terrain_materials` | `Material[]` | `[]` | — | Array of Godot `Material` resources. Index = material ID assigned by `UnderGenMaterialStamperNode`. |
| `liquid_material` | `Material` | `null` | — | Material for liquid/water surfaces. |
| `compute_shader` | `RDShaderFile` | `null` | — | Compute shader for optimized mesh generation (optional). |

---

### 3.13 UnderGenSceneSpawnerNode

**Editor Label**: Scene Spawner
**Role**: **Terminal node.** Instantiates a **single** `PackedScene` at every point in the input `PointSet`. Requires `UnderGenWorld3D` as parent.

> **Per-zone spawning**: To spawn different scenes per zone, place `UnderGenPointFilterNode` instances upstream — each filtering for a specific zone — and wire each to its own spawner. Use `zone_match_mode: 2` (Exclude) on a final filter to catch all remaining zones.

| Port | Direction | Type | Description |
|:-----|:----------|:-----|:------------|
| 0 | **in** | `3` Point Set | Points to spawn at. |
| 0 | **out** | `3` Point Set | Pass-through for chaining. |

| Property | Type | Default | Range | Description |
|:---------|:-----|:--------|:------|:------------|
| `scene_to_spawn` | `PackedScene` | `null` | — | The scene to instantiate at each point. Accepts `res://` path strings in JSON. |
| `spawn_probability` | `float` | `1.0` | `0..1` | Probability each point spawns. Also multiplied by point's `density`. |
| `random_y_rotation` | `bool` | `true` | — | Randomize Y-axis rotation of spawned instances. |
| `random_seed` | `int` | `0` | — | Seed for RNG. |

---

## 4. Full Examples

### 4.1 Complete Grammar JSON

```json
{
  "axiom": "Start",
  "state_variables": ["keys", "danger_level"],
  "room_types": [
    { "symbol": "Entry",     "color": [0.0, 1.0, 1.0, 1.0], "min_size": [4, 3, 4], "max_size": [8, 5, 8] },
    { "symbol": "Hallway",   "color": [1.0, 0.6, 0.0, 1.0], "min_size": [3, 3, 6], "max_size": [6, 5, 12] },
    { "symbol": "BossRoom",  "color": [1.0, 0.0, 0.0, 1.0], "min_size": [8, 5, 8], "max_size": [14, 8, 14], "vox_path": "res://vox/boss_arena.vox" },
    { "symbol": "Treasure",  "color": [1.0, 0.84, 0.0, 1.0], "min_size": [3, 3, 3], "max_size": [6, 5, 6] },
    { "symbol": "KeyRoom",   "color": [0.0, 0.8, 0.4, 1.0], "min_size": [4, 3, 4], "max_size": [7, 5, 7] },
    { "symbol": "MobArena",  "color": [0.8, 0.2, 0.8, 1.0], "min_size": [5, 4, 5], "max_size": [10, 7, 10] }
  ],
  "rules": [
    {
      "rule_name": "Root Split",
      "lhs_symbol": "Start",
      "probability": 1.0,
      "entry_node_id": "entry",
      "exit_node_id": "boss",
      "rhs_nodes": [
        { "id": "entry", "symbol": "Entry",   "constraints": { "fixed_pos": [0, 0, -20] } },
        { "id": "hall",  "symbol": "Hallway" },
        { "id": "boss",  "symbol": "BossRoom", "constraints": { "fixed_pos": [0, 0, 20] } }
      ],
      "rhs_edges": [
        { "from": "entry", "to": "hall", "type": "corridor" },
        { "from": "hall",  "to": "boss", "type": "corridor" }
      ]
    },
    {
      "rule_name": "Hallway → Key + Treasure",
      "lhs_symbol": "Hallway",
      "probability": 1.0,
      "entry_node_id": "key",
      "exit_node_id": "treasure",
      "rhs_nodes": [
        { "id": "key",      "symbol": "KeyRoom" },
        { "id": "treasure", "symbol": "Treasure" }
      ],
      "rhs_edges": [
        { "from": "key", "to": "treasure", "type": "corridor" }
      ],
      "actions": [{ "var": "keys", "delta": 1 }]
    },
    {
      "rule_name": "Hallway → Mob Arena (conditional)",
      "lhs_symbol": "Hallway",
      "probability": 0.5,
      "condition_var": "danger_level",
      "condition_op": ">=",
      "condition_val": 2,
      "rhs_nodes": [
        { "id": "arena", "symbol": "MobArena" }
      ],
      "rhs_edges": [],
      "actions": [{ "var": "danger_level", "delta": 1 }]
    }
  ]
}
```

### 4.2 Complete Pipeline JSON

```json
{
  "nodes": [
    {
      "name": "Grammar",
      "type": "UnderGenGrammarNode",
      "pos": [100, 200],
      "properties": {
        "grammar_resource_path": "res://grammars/dungeon.tres",
        "iterations": 5,
        "max_nodes": 150
      }
    },
    {
      "name": "BSP_Placer",
      "type": "UnderGenBSPPlacerNode",
      "pos": [350, 200],
      "properties": {
        "grid_size_x": 64,
        "grid_size_y": 20,
        "grid_size_z": 64
      }
    },
    {
      "name": "Noise_Blend",
      "type": "UnderGenNoiseNode",
      "pos": [500, 150],
      "properties": {
        "noise_scale": 40.0,
        "noise_intensity": 0.35,
        "noise_frequency": 0.025,
        "noise_seed": 42
      }
    },
    {
      "name": "Bezier_Carver",
      "type": "UnderGenBezierCarverNode",
      "pos": [600, 200],
      "properties": {
        "path_brush_min_radius": 2,
        "path_brush_max_radius": 5,
        "path_segments": 2,
        "path_bend_factor": 0.5,
        "path_wobble_magnitude": 0.4,
        "cave_ruggedness": 1.2,
        "cave_floor_ruggedness": 0.3,
        "floor_flattening": 0.6
      }
    },
    {
      "name": "Smooth",
      "type": "UnderGenSmoothNode",
      "pos": [750, 200],
      "properties": {
        "smoothing_strength": 1
      }
    },
    {
      "name": "Material_Stamp",
      "type": "UnderGenMaterialStamperNode",
      "pos": [850, 200],
      "properties": {
        "default_material_id": 0,
        "zone_material_map": [
          { "zone_name": "bossroom", "material_id": 1 },
          { "zone_name": "corridor", "material_id": 0 },
          { "zone_name": "treasure", "material_id": 2 }
        ]
      }
    },
    {
      "name": "Surface_Sampler",
      "type": "UnderGenSurfaceSamplerNode",
      "pos": [950, 150],
      "properties": {
        "surface_type": 0,
        "zone_filter": "mobarena,bossroom"
      }
    },
    {
      "name": "Spawn_Filter_Boss",
      "type": "UnderGenPointFilterNode",
      "pos": [1050, 140],
      "properties": {
        "required_zone_name": "BossRoom",
        "zone_match_mode": 0,
        "min_spacing": 4.0,
        "min_slope": 0.0,
        "max_slope": 0.3
      }
    },
    {
      "name": "Spawner_Boss",
      "type": "UnderGenSceneSpawnerNode",
      "pos": [1180, 140],
      "properties": {
        "scene_to_spawn": "res://entities/boss.tscn",
        "spawn_probability": 1.0,
        "random_seed": 100
      }
    },
    {
      "name": "Spawn_Filter_Default",
      "type": "UnderGenPointFilterNode",
      "pos": [1050, 200],
      "properties": {
        "required_zone_name": "BossRoom",
        "zone_match_mode": 2,
        "min_spacing": 3.0,
        "min_slope": 0.0,
        "max_slope": 0.3
      }
    },
    {
      "name": "Spawner_Default",
      "type": "UnderGenSceneSpawnerNode",
      "pos": [1180, 200],
      "properties": {
        "scene_to_spawn": "res://entities/enemy.tscn",
        "spawn_probability": 0.8,
        "random_seed": 200
      }
    },
    {
      "name": "Mesher",
      "type": "UnderGenMesherNode",
      "pos": [950, 300],
      "properties": {
        "chunk_size": 16,
        "voxel_size": 1.0,
        "generate_collision": true
      }
    }
  ],
  "connections": [
    { "from_node": "Grammar",          "from_port": 0, "to_node": "BSP_Placer",        "to_port": 1 },
    { "from_node": "BSP_Placer",       "from_port": 0, "to_node": "Noise_Blend",       "to_port": 0 },
    { "from_node": "Noise_Blend",      "from_port": 0, "to_node": "Bezier_Carver",     "to_port": 0 },
    { "from_node": "Bezier_Carver",    "from_port": 0, "to_node": "Smooth",            "to_port": 0 },
    { "from_node": "Smooth",           "from_port": 0, "to_node": "Material_Stamp",    "to_port": 0 },
    { "from_node": "Material_Stamp",   "from_port": 0, "to_node": "Surface_Sampler",   "to_port": 0 },
    { "from_node": "Material_Stamp",   "from_port": 0, "to_node": "Mesher",            "to_port": 0 },
    { "from_node": "Surface_Sampler",       "from_port": 1, "to_node": "Spawn_Filter_Boss",   "to_port": 0 },
    { "from_node": "Spawn_Filter_Boss",     "from_port": 0, "to_node": "Spawner_Boss",         "to_port": 0 },
    { "from_node": "Surface_Sampler",       "from_port": 1, "to_node": "Spawn_Filter_Default", "to_port": 0 },
    { "from_node": "Spawn_Filter_Default",  "from_port": 0, "to_node": "Spawner_Default",      "to_port": 0 }
  ]
}
```

### 4.3 Minimal Pipeline (Grammar → Placer → Carver → Mesher)

```json
{
  "nodes": [
    {
      "name": "Grammar",
      "type": "UnderGenGrammarNode",
      "properties": { "grammar_resource_path": "res://grammars/simple.tres", "iterations": 3 }
    },
    {
      "name": "Placer",
      "type": "UnderGenBSPPlacerNode",
      "properties": { "grid_size_x": 48, "grid_size_y": 16, "grid_size_z": 48 }
    },
    {
      "name": "Carver",
      "type": "UnderGenAStarCarverNode",
      "properties": { "path_brush_min_radius": 2, "path_brush_max_radius": 4 }
    },
    {
      "name": "Mesher",
      "type": "UnderGenMesherNode",
      "properties": { "chunk_size": 16 }
    }
  ],
  "connections": [
    { "from_node": "Grammar", "from_port": 0, "to_node": "Placer", "to_port": 1 },
    { "from_node": "Placer",  "from_port": 0, "to_node": "Carver", "to_port": 0 },
    { "from_node": "Carver",  "from_port": 0, "to_node": "Mesher", "to_port": 0 }
  ]
}
```

---

## 5. API Usage

### GDScript

```gdscript
# ── Load grammar from JSON ──────────────────────────────────
var grammar = LevelGrammarSpec.load_from_json_file("res://grammars/dungeon.grammar.json")
grammar.save_to_file("res://grammars/dungeon.tres")   # also works with editor

# ── Load pipeline from JSON ─────────────────────────────────
var pipeline = UnderGenPipeline.load_from_json_file("res://pipelines/dungeon.pipeline.json")
pipeline.save_to_file("res://pipelines/dungeon.tres")  # also works with editor

# ── Execute ──────────────────────────────────────────────────
var world = UnderGenWorld3D.new()
add_child(world)
world.set_pipeline(pipeline)

var result = {}
var ok = pipeline.execute_pipeline({"seed": 12345}, result)
if ok:
    print("Generation succeeded! Result keys: ", result.keys())
```

### C++

```cpp
// Load grammar from JSON
Ref<LevelGrammarSpec> grammar = LevelGrammarSpec::load_from_json_file("res://grammars/dungeon.grammar.json");

// Load pipeline from JSON
Ref<UnderGenPipeline> pipeline = UnderGenPipeline::load_from_json_file("res://pipelines/dungeon.pipeline.json");

// Execute
Dictionary result;
bool ok = pipeline->execute_pipeline(pair("seed", 12345), result);
```

### JSON → .tres conversion (for editor compatibility)

```gdscript
# Batch convert all .grammar.json files to .tres
static func convert_grammars():
    var dir = DirAccess.open("res://grammars/")
    dir.list_dir_begin()
    var file = dir.get_next()
    while file != "":
        if file.ends_with(".grammar.json"):
            var spec = LevelGrammarSpec.load_from_json_file("res://grammars/" + file)
            var tres_path = "res://grammars/" + file.replace(".grammar.json", ".tres")
            spec.save_to_file(tres_path)
            print("Converted: ", file, " → ", tres_path.get_file())
        file = dir.get_next()
```

---

## Type Quick Reference

| JSON type | Godot type | Example |
|:----------|:-----------|:--------|
| `string` | `String` | `"hello"` |
| `float` / `int` | `float` / `int` | `3.14`, `42` |
| `bool` | `bool` | `true`, `false` |
| `[x, y]` | `Vector2` | `[100, 200]` |
| `[x, y, z]` (int) | `Vector3i` | `[5, 3, 5]` |
| `[x, y, z]` (float) | `Vector3` | `[1.5, 0, -2.0]` |
| `[r, g, b, a]` | `Color` | `[1, 0, 0, 1]` (red) |
| `[a, b, c]` | `Array` | `[1, 2, 3]` |
| `{ "k": v }` | `Dictionary` | `{ "key": 1 }` |
| `"res://..."` | `Resource` path | `"res://scenes/enemy.tscn"` |
| `null` | null / empty | skip the property or set `null` |

---

*Document generated from the UnderGen GDExtension source code. For questions or contributions, refer to the project repository.*
