# Spatial generation architecture

UnderGen 2.0 separates level generation into three authored stages:

1. **Topology** produces an `UnderGenSemanticGraph`. Nodes express spatial roles and gameplay tags independently; edges express traversal, visibility, enclosure, exposure, wrapping, crossing, containment, and vertical relationships.
2. **Embedding** produces an `UnderGenEmbeddedLayout`. It contains mutable 3D spaces, stateful paths, coarse semantic fields, constraint results, a revision, and dirty regions.
3. **Geometry** first produces an ordered `UnderGenGeometryPlan`, then realizes that plan into the existing `DensityGrid`. The realized generation context retains all three staged artifacts and also emits legacy `rooms` and `edges` arrays for existing downstream nodes.

The recommended pipeline is:

```text
Topology Builder
  -> Spatial Embedder
  -> Spatial Validator
  -> Geometry Planner
  -> Geometry Realizer
  -> existing detail/material/sampling nodes
  -> Mesher
```

The pipeline editor's **Spatial Template** menu creates this stack using Layered Chasm, Compact Cave, or Branching Cavern. Presets remain ordinary `UnderGenGenerationPreset` resources and are editable in the Godot Inspector.

## Building a playable spatial level

There are two spatial authoring modes:

- The composition templates generate terrain and semantic spaces directly. They are useful for shaping a cavern, but do not invent project-specific player, portal, or encounter conventions.
- **Playable from Grammar (Layered Chasm)** is the compatibility workflow. It keeps an existing `LevelGrammarResource` as the gameplay topology, embeds that graph in 3D, and then realizes and decorates it through the spatial stack.

The playable hybrid pipeline is:

```text
Grammar Expander
  -> Topology Builder (legacy input)
  -> Spatial Embedder
  -> Spatial Validator
  -> Geometry Planner
  -> Geometry Realizer
  -> geological/detail passes
  -> VOX Stamper
  -> Gameplay Markers
  -> sampler/spawner passes
  -> Mesher
```

`UnderGenTopologyBuilderNode` preserves each grammar node's `vox_path`, constraints, and smoothing/warping exclusions as topology parameters. The embedder carries those parameters into its spaces, and the geometry realizer publishes them again in the legacy `rooms` array. Existing VOX, material, detail, sampler, and spawner nodes can therefore remain downstream of the new architecture.

`UnderGenGameplayMarkerNode` is the bridge to `LevelGenerator`. It emits one `PlayerStart` from the semantic entry space and `exit_portal` records from semantic exit spaces. By default it replaces the equivalent route markers found inside the entry/exit VOX models, while preserving unrelated authored markers such as chests. This prevents duplicate player or portal spawns and makes those route markers follow live embedding edits.

To use it in the editor:

1. Open the pipeline editor and choose **Spatial Template -> Playable from Grammar (Layered Chasm)**.
2. Choose the grammar resource when prompted. Existing grammars such as `small_grammar.tres` work without conversion.
3. Add or copy the biome-specific material, detail, sampler, and scene-spawner passes used by the old pipeline.
4. Save the new pipeline and assign it to the biome or `UnderGenWorld3D` in the same way as an old pipeline. Keep the grammar assigned on the world/biome if runtime grammar overrides are required.
5. Generate once, select `UnderGenWorld3D`, and use the viewport handles to adjust the retained embedding. Player and exit records are refreshed when spatial geometry is rebuilt.

`res://resources/undergen/pipelines/forest_playable_spatial.tres` is a complete example based on `small_grammar.tres` and the forest pipeline's VOX palette, materials, tree pass, and wolf pass. It is intentionally a separate resource, so existing playable pipelines keep their current generation behavior until explicitly migrated.

## Runtime embedding edits

The latest context is retained by `UnderGenWorld3D`. Runtime code can edit one space or a complete elevation band:

```gdscript
# Move the authored middle traversal band six voxels lower.
$UnderGenWorld3D.move_elevation_band(0, -6.0)

# Move one embedded landmark to an exact position.
$UnderGenWorld3D.move_embedded_space("upper_overlook", Vector3(62, 54, 38))

# Separate the edit from realization when batching several changes.
$UnderGenWorld3D.set_embedded_space_elevation("upper_overlook", 48.0, false)
$UnderGenWorld3D.rebuild_spatial_geometry()
```

Moving a space retargets connected path endpoints. Moving a band translates internal path controls when both endpoints belong to the band and retargets only the affected endpoint otherwise. The layout records the union of old and new bounds. Geometry is replanned, ordered CSG is replayed only in the dirty volume, and only intersecting terrain chunks are remeshed.

Signals `spatial_layout_changed` and `spatial_geometry_rebuilt` let gameplay/detail systems invalidate their own derived data. Detail passes that are not represented in the geometry plan should use those signals to restamp affected content.

### Editor handles

After generating a spatial pipeline in the editor, select the `UnderGenWorld3D` node to edit its retained embedding directly in the 3D viewport. Selecting or remeshing the world stays in the 3D editor; the UnderGen main screen opens only for grammar and pipeline resources.

- Use the 3D viewport's **Transform Mode**, not the separate selection-only mode introduced in Godot 4.6.
- Drag a cyan cross at a box center to move one embedded space in the camera-facing plane. The viewport uses a 22-pixel screen-space hit area, even when Godot does not render the optional ring icon. The wireframe box itself is visualization only.
- Shift-drag at the center of an orange vertical band marker to move that entire elevation band along embedding Y. The orange path skeleton itself is visualization only.
- Releasing a handle rebuilds the accumulated dirty region once. Pressing Escape cancels the drag, and the completed operation participates in editor Undo/Redo.

The viewport drag deliberately previews only the embedding boxes and paths. Density replay and chunk remeshing happen on release, keeping large generated levels responsive while a handle is moving.

## Composition controls

The spatial validator reports:

- required traversal and endpoint validity;
- containment and vertical-order constraints;
- positive-mass to primary-void ratio;
- boundary-route and cross-void route ratios;
- visible elevation-band count.

Embedding also builds eight low-resolution fields: openness, verticality, exposure, enclosure, occlusion, prominence, connectivity pressure, and surface suitability. Geometry planning samples these fields when sizing path-state segments. They remain independent of the density grid so other stages can use them without reading voxels.

Paths carry ordered state keys rather than a single tunnel profile. Width, height, floor flatness, lateral bias, wall presence, exposure, and local noise scale are interpolated into separate sweep operations. Boundary routes add ledges and optional undercuts; crossings add a positive bridge span; all operations preserve source IDs for later rebuilding.

## Compatibility and trade-offs

Existing grammar/BSP and outdoor pipelines remain valid. Grammar output can be adapted through `UnderGenTopologyBuilderNode`, and `UnderGenGeometryRealizerNode` emits the legacy context expected by current detail, material, sampler, spawner, and mesher nodes. Seed input now reaches disconnected source nodes through explicit pipeline defaults.

The staged workflow adds concepts, but reduces per-level graph work. A level designer normally chooses a preset, edits semantic roles/relations and a handful of composition values, then inspects validator output and the 3D gizmo. Fine voxel operations remain reusable generator infrastructure instead of being authored per level.

The refactor intentionally stops geometry from deciding topology. Direct voxel edits can still be performed downstream, but they cannot create or repair semantic connectivity automatically. Likewise, runtime dirty replay reconstructs the planned macro/route geometry; downstream decoration systems remain responsible for refreshing their own derived objects in the signaled dirty regions.
