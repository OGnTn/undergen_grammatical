class_name WorldSpawnableObject extends Resource

@export_group("Scene Data")
@export var scene: PackedScene

@export_group("Placement Rules")
@export var align_to_surface_normal: bool = false
@export var random_y_rotation: bool = true
@export_range(0.0, 1.0) var spawn_chance: float = 0.5
@export var scan_bounds: AABB = AABB(Vector3(-0.5, -0.5, -0.5), Vector3(1, 1, 1))
@export var position_offset: Vector3 = Vector3.ZERO

@export_group("Slope Constraints")
@export_range(-1.0, 1.0) var normal_y_min: float = 0.8
@export_range(-1.0, 1.0) var normal_y_max: float = 1.0

@export_group("Zone Rules (Director)")
## The exact Zone Name required to spawn here (e.g., 'safe_glade', 'spider_tunnel').
## Matches the 'type' string defined in your Grammar.
## Leave empty to spawn in any zone.
@export var required_zones: Array[String] = [""]

## If true, this object will only spawn if the Zone is NOT a Room (i.e., it is a Corridor/Path).
## Note: This relies on your Grammar correctly tagging paths vs rooms.
@export var path_only: bool = false

@export_group("Physics")
@export var check_collision: bool = true
