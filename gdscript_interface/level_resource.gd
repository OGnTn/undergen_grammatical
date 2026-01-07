class_name LevelResource extends Resource

@export_group("Logic")
@export var grammar: LevelGrammarResource

@export_group("Physics (Global Defaults)")
@export var world_size_x: int = 1024
@export var world_size_y: int = 64
@export var world_size_z: int = 1024
@export var chunk_size: int = 16
@export var cube_size: int = 1
## Maximum number of logical nodes allowed in the grammar. Stops expansion if reached.
@export var max_rooms: int = 100
@export var room_size_range_x: Vector2 = Vector2(10, 20)
@export var room_size_range_y: Vector2 = Vector2(5, 10)
@export var room_size_range_z: Vector2 = Vector2(10, 20)
@export var dungeon_mode: bool = false
@export var path_connect_from_ground_level: bool = true

@export_subgroup("Bezier Paths")
@export_range(1, 10) var path_segments: int = 5
@export_range(0.0, 2.0) var path_bend_factor: float = 0.5
@export_range(0.0, 10.0) var path_wobble_magnitude: float = 0.5
@export_range(0.0, 1.0) var path_wobble_frequency: float = 0.1

@export_subgroup("Brush Settings")
@export_range(1, 10) var path_brush_min_radius: int = 2
@export_range(1, 10) var path_brush_max_radius: int = 4
@export var path_use_square_brush: bool = false

@export_subgroup("Smoothing")
@export var smoothing_enabled: bool = false
@export_range(1, 5) var smoothing_strength: int = 1

@export_group("Art & Materials")
## Material ID 0 (Default)
@export var base_material: Material
## Material IDs 1+ (Mapped by index)
@export var custom_materials: Array[Material]
@export var liquid_material: Material

@export_group("Gameplay (Spawning)")
@export var spawn_registry: Array[WorldSpawnableObject]

@export_group("Liquid Simulation")
@export var liquid_source_count: int = 10
@export var liquid_spread_range: int = 15

@export_group("Vox Integration")
## Maps Palette Index (int 0-255) to a WorldSpawnableObject
@export var vox_spawn_map: Dictionary = {}
## Maps Palette Index (int 0-255) to a Material (overrides generated material)
@export var vox_material_map: Dictionary = {}

@export_group("Legacy / Internal")
@export var noise_texture: FastNoiseLite
@export var decorations: Dictionary[PackedScene, float] = {}
@export var entities: Dictionary[PackedScene, float] = {}
@export var light_color: Color = Color(1, 1, 1, 1)
@export var light_energy: float = 1.0
