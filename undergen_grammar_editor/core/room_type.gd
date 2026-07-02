@tool
class_name RoomType extends Resource
## A single room type entry in the grammar palette.
## Defines the visual and spatial identity of one symbol.

@export var symbol: String = "Room"
@export var color: Color = Color(0.4, 0.6, 0.9)
@export var weight: float = 1.0

@export_group("Spatial Size")
@export var min_size: Vector3i = Vector3i(5, 3, 5)
@export var max_size: Vector3i = Vector3i(10, 6, 10)

@export_group("Asset")
@export var vox_path: String = ""
@export var exclude_from_smoothing: bool = false
