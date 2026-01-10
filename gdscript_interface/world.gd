class_name World extends Node3D


var world_size_x: int
var world_size_y: int
var world_size_z: int

var chunk_size: int
var cube_size: float

var chunk_count_x: int
var chunk_count_y: int
var chunk_count_z: int

var chunks: Array[MCChunk] = []
var liquid_chunks: Array[MCChunk] = []
var surface_normals: Dictionary = {}
var face_normals: Dictionary[Vector3, Vector3] = {}

@export_group("Data")
@export var SOLID: float = 0.0
@export var OPEN: float = 1.0
@export var surface_threshold: float = 0.05

@export_group("Generation")
@export var level_resource: LevelResource
@export var density_grid_resource: LevelDensityGrid

@export_group("Scene Nodes")
@export var poisson_disk_sampler: PoissonDiskSampler
@export var generator: Generator
@export var populator: Node3D
@export var water_sim: WaterSimulator
@export var navigation_region: NavigationRegion3D

@export_category("Debug")
@export var generate: bool:
	set(value):
		generate_world()

@export var debug_show_zones: bool = false:
	set(value):
		debug_show_zones = value
		_update_zone_visualizer()

var _zone_visualizer = null

func get_closest_ground_position(position: Vector3) -> Vector3:
	return density_grid_resource.get_closest_ground_position(position)

func _update_zone_visualizer():
	if debug_show_zones:
		if not _zone_visualizer:
			var script = load("res://systems/debug/zone_visualizer.gd")
			if script:
				_zone_visualizer = Node3D.new()
				_zone_visualizer.set_script(script)
				_zone_visualizer.name = "ZoneVisualizer"
				add_child(_zone_visualizer)
				_zone_visualizer.world = self
		_zone_visualizer.visualize()
	else:
		if _zone_visualizer:
			_zone_visualizer.queue_free()
			_zone_visualizer = null

func export_scene():
	generate_world()

func init_level_resource():
	world_size_x = level_resource.world_size_x
	world_size_y = level_resource.world_size_y
	world_size_z = level_resource.world_size_z
	chunk_size = level_resource.chunk_size
	cube_size = level_resource.cube_size
	chunk_count_x = ceil(float(world_size_x) / chunk_size)
	chunk_count_y = ceil(float(world_size_y) / chunk_size)
	chunk_count_z = ceil(float(world_size_z) / chunk_size)
	
	if not density_grid_resource:
		printerr("ERROR: LevelDensityGrid resource is not assigned to World.density_grid_resource!")
		return
	
	density_grid_resource.initialize_grid(world_size_x, world_size_y, world_size_z, SOLID)
	density_grid_resource.path_dungeon_mode = level_resource.dungeon_mode
	density_grid_resource.path_connect_from_ground_level = level_resource.path_connect_from_ground_level
	density_grid_resource.path_segments = level_resource.path_segments
	density_grid_resource.path_bend_factor = level_resource.path_bend_factor
	density_grid_resource.path_wobble_magnitude = level_resource.path_wobble_magnitude
	density_grid_resource.path_wobble_frequency = level_resource.path_wobble_frequency
	
	density_grid_resource.path_brush_min_radius = level_resource.path_brush_min_radius
	density_grid_resource.path_brush_max_radius = level_resource.path_brush_max_radius
	density_grid_resource.path_use_square_brush = level_resource.path_use_square_brush
	density_grid_resource.smoothing_enabled = level_resource.smoothing_enabled
	density_grid_resource.smoothing_strength = level_resource.smoothing_strength
	
	density_grid_resource.surface_threshold = surface_threshold
	density_grid_resource.room_min_size = Vector3i(level_resource.room_size_range_x[0], level_resource.room_size_range_y[0], level_resource.room_size_range_z[0])
	density_grid_resource.room_max_size = Vector3i(level_resource.room_size_range_x[1], level_resource.room_size_range_y[1], level_resource.room_size_range_z[1])

func generate_world() -> void:
	populator.finished_spawning.connect(_on_all_objects_spawned)
	print(str(multiplayer.get_unique_id()) + ": " + "Generating level")
	init_level_resource()
	init_chunks()
	print(str(multiplayer.get_unique_id()) + ": " + "initialized chunks")
	# 1. Generate Terrain Density
	# 3. Generate Terrain (The Integration Point)
	# Check if we have a valid Graph Resource assigned in the Inspector
	# Check for Grammar in LevelResource
	if level_resource.grammar:
		print(str(multiplayer.get_unique_id()) + ": " + "Generating from Procedural Grammar (LevelResource)...")
		
		# 2. Run Grammar
		var gen = GrammarGenerator.new(level_resource.grammar, 1337)
		var result = gen.generate(4, {}, level_resource.max_rooms)
		
		# 3. Pass to C++
		density_grid_resource.generate_from_graph(
			result["nodes"],
			result["edges"],
			cube_size,
			1337,
			{
				"vox_spawn_map": level_resource.vox_spawn_map,
				"vox_material_map": level_resource.vox_material_map
			}
		)
	else:
		# Fallback: Use the legacy random generation if no grammar is provided
		print(str(multiplayer.get_unique_id()) + ": " + "Generating Default Noise Level...")
		generator.generate_level(10, density_grid_resource)
	print(str(multiplayer.get_unique_id()) + ": " + "Generated level")
	surface_normals = density_grid_resource.surface_normals
	
	# Run simulation (C++) and get the new DensityGrid specifically for water
	print(str(multiplayer.get_unique_id()) + ": " + "Generating liquid")
	var liquid_grid_resource = density_grid_resource.generate_liquid_grid()
	print(str(multiplayer.get_unique_id()) + ": " + "Generated liquid")
	liquid_grid_resource.surface_threshold = surface_threshold
	
	# Assign the liquid grid to all liquid chunks
	for l_chunk in liquid_chunks:
		l_chunk.set_density_grid(liquid_grid_resource)
	
	print(str(multiplayer.get_unique_id()) + ": " + "building mesh")
	build_mesh()
	print(str(multiplayer.get_unique_id()) + ": " + "built mesh")
	if (multiplayer.is_server()):
		navigation_region.bake_navigation_mesh()
	populator.prepare_spawn_list(level_resource.spawn_registry)
	
	_update_zone_visualizer()

func _on_all_objects_spawned():
	bake_gi()
	Network.on_player_loaded_level.rpc_id(1)

func bake_gi():
	#await get_tree().create_timer(1.5).timeout
	print(str(multiplayer.get_unique_id()) + ": " + "baking GI")
	%VoxelGI.position = Vector3(world_size_x * cube_size / 2, world_size_y * cube_size / 2, world_size_z * cube_size / 2)
	%VoxelGI.size = Vector3(world_size_x * cube_size, world_size_y * cube_size, world_size_z * cube_size)
	%VoxelGI.bake()


func get_true_surface_position(approx_pos_ws: Vector3) -> Vector3:
	# 1. Find which chunk this world-space position is in.
	var chunk_x = floori(approx_pos_ws.x / (chunk_size * cube_size))
	var chunk_y = floori(approx_pos_ws.y / (chunk_size * cube_size))
	var chunk_z = floori(approx_pos_ws.z / (chunk_size * cube_size))
	
	var chunk_index = to_index(chunk_x, chunk_y, chunk_z)
	if chunk_index < 0 or chunk_index >= chunks.size():
		return approx_pos_ws # Return approximation if out of bounds
		
	var chunk: MCChunk = chunks[chunk_index]
	if not chunk or not is_instance_valid(chunk.mesh):
		return approx_pos_ws

	# 2. Get the mesh vertex data.
	var mesh: ArrayMesh = chunk.mesh
	if mesh.get_surface_count() == 0:
		return approx_pos_ws
	var surface_arrays = mesh.surface_get_arrays(0)
	var vertices: PackedVector3Array = surface_arrays[Mesh.ARRAY_VERTEX]
	
	if vertices.is_empty():
		return approx_pos_ws

	# 3. Find the closest vertex on that mesh.
	var closest_vertex_ws: Vector3
	var min_dist_sq = -1.0
	
	# The chunk's vertices are in its local space, so we must add the chunk's
	# own position to get them into world space for comparison.
	var chunk_pos_ws = chunk.global_position

	for v in vertices:
		var vertex_ws = v + chunk_pos_ws
		var dist_sq = approx_pos_ws.distance_squared_to(vertex_ws)
		
		if min_dist_sq < 0 or dist_sq < min_dist_sq:
			min_dist_sq = dist_sq
			closest_vertex_ws = vertex_ws
			
	return closest_vertex_ws

func align_with_y(xform, new_y):
	var up_direction = Vector3.UP
	var rotation_quat = Quaternion(up_direction, new_y)
	xform.basis = Basis(rotation_quat).orthonormalized()
	return xform

func get_n_closest_surface_points(x, y, z, n):
	var points = surface_normals.keys()
	points.sort_custom(func(a, b): a.distance_to(Vector3i(x, y, z)) < b.distance_to(Vector3i(x, y, z)))
	return points.slice(0, n)

func get_closest_surface_point(x, y, z):
	var closest_point: Vector3i
	var closest_dist: float
	for p in surface_normals.keys():
		if (!closest_point):
			closest_point = p
			closest_dist = p.distance_to(Vector3i(x, y, z))
		else:
			var dist = p.distance_to(Vector3i(x, y, z))
			if dist < closest_dist:
				closest_point = p
				closest_dist = dist
	return closest_point

func mark_brush(x, y, z, radius_lower, radius_higher, val):
	density_grid_resource._mark_brush(Vector3i(x, y, z), radius_lower, radius_higher, val)

func modify_terrain_material(position: Vector3, radius: int, material_index: int):
	# 1. Convert World Position to Grid Coordinates
	var cx = int(round(position.x / cube_size))
	var cy = int(round(position.y / cube_size))
	var cz = int(round(position.z / cube_size))
	
	var r_sq = float(radius * radius)
	var center = Vector3(cx, cy, cz)
	
	# 2. Iterate through the brush area
	for x in range(cx - radius, cx + radius + 1):
		for y in range(cy - radius, cy + radius + 1):
			for z in range(cz - radius, cz + radius + 1):
				# Check distance for sphere shape
				var current_pos = Vector3(x, y, z)
				if current_pos.distance_squared_to(center) <= r_sq:
					var grid_pos = Vector3i(x, y, z)
					
					# 3. Apply the Material ID
					# This calls the C++ method we added to DensityGrid earlier
					density_grid_resource.set_material_id(grid_pos, material_index)

	# 4. Identify and Update Affected Chunks
	# We can reuse the existing helper for this
	_update_chunks_in_bounds(cx, cy, cz, radius)

func init_chunks() -> void:
	# 1. Clear Old Chunks
	# Iterating get_children() in GDScript can also be slow if there are 4000+ nodes.
	# It's better to clear the arrays we know about.
	for c in chunks:
		if is_instance_valid(c): c.queue_free()
	for c in liquid_chunks:
		if is_instance_valid(c): c.queue_free()
		
	chunks.clear()
	liquid_chunks.clear()
	
	# 2. Initialize using C++ Manager
	var chunk_manager = ChunkManager.new()
	#var compute_shader = load("res://assets/mc_compute_multimat.glsl")
	var compute_shader = null
	
	print(str(multiplayer.get_unique_id()) + ": " + "Batch init chunks (C++)")

	# 1. Prepare Materials Array
	var materials = [level_resource.base_material]
	if level_resource.custom_materials:
		materials.append_array(level_resource.custom_materials)
		
	# This single call replaces the entire nested loop
	chunk_manager.init_chunks_batch(
		self, # Parent node (World)
		chunks, # Array to fill (Terrain)
		liquid_chunks, # Array to fill (Liquid)
		density_grid_resource,
		materials, # Pass array of materials
		level_resource.liquid_material,
		compute_shader,
		chunk_size,
		cube_size,
		chunk_count_x,
		chunk_count_y,
		chunk_count_z
	)
	
	chunk_manager.free() # We don't need the manager node anymore

func build_mesh() -> void:
	print(str(multiplayer.get_unique_id()) + ": " + "Building mesh (C++ MCChunk)")
	var chunk_count = 0
	
	# Build Terrain Meshes
	for chunk in chunks:
		chunk.generate_mesh_from_density_grid()
		chunk_count += 1
		
	# Build Liquid Meshes [NEW]
	for l_chunk in liquid_chunks:
		l_chunk.generate_mesh_from_density_grid()
		
	print(str(multiplayer.get_unique_id()) + ": " + "Built ", chunk_count, " terrain chunks + liquid.")

func is_in_bounds(x, y, z):
	return x > 1 and x < world_size_x - 1 and y > 1 and y < world_size_y - 1 and z > 1 and z < world_size_z - 1

func get_cell(x: int, y: int, z: int) -> float:
	return density_grid_resource.get_cell(Vector3i(x, y, z), SOLID)

func set_cell(x: int, y: int, z: int, val: float):
	density_grid_resource.set_cell(Vector3i(x, y, z), val)
	var chunk_x = floor(float(x) / chunk_size)
	var chunk_y = floor(float(y) / chunk_size)
	var chunk_z = floor(float(z) / chunk_size)
	
	var affected_chunk: MCChunk = chunks[to_index(chunk_x, chunk_y, chunk_z)]
	if affected_chunk:
		pass

func to_chunk_index(x: int, y: int, z: int) -> int:
	return x + y * chunk_size + z * chunk_size * chunk_size

func to_index(x: int, y: int, z: int) -> int:
	return x + y * chunk_count_x + z * chunk_count_x * chunk_count_y


func modify_terrain_additive(position: Vector3, radius: int, is_digging: bool):
	# 1. Convert World Position to Grid Coordinates
	var cx = int(round(position.x / cube_size))
	var cy = int(round(position.y / cube_size))
	var cz = int(round(position.z / cube_size))
	
	var r_sq = float(radius * radius)
	var center = Vector3(cx, cy, cz)
	
	# 2. Iterate through the brush area
	for x in range(cx - radius, cx + radius + 1):
		for y in range(cy - radius, cy + radius + 1):
			for z in range(cz - radius, cz + radius + 1):
				# Check distance for sphere shape
				var current_pos = Vector3(x, y, z)
				if current_pos.distance_squared_to(center) <= r_sq:
					var grid_pos = Vector3i(x, y, z)
					
					# Get current value (defaulting to SOLID/0.0 if uninitialized, 
					# but usually you want to default to OPEN/1.0 for empty space)
					# Assuming default is OPEN based on your context:
					var current_val = density_grid_resource.get_cell(grid_pos, OPEN)
					
					# 3. Add or Subtract
					# If digging, we want to go towards 1.0 (Add)
					# If building, we want to go towards 0.0 (Subtract)
					if is_digging:
						current_val += .01
					else:
						current_val -= .01
					
					# 4. Clamp to valid range [0.0, 1.0]
					current_val = clamp(current_val, SOLID, OPEN)
					
					density_grid_resource.set_cell(grid_pos, current_val)
	
	# 5. Identify and Update Affected Chunks (Standard Logic)
	_update_chunks_in_bounds(cx, cy, cz, radius)

func _update_chunks_in_bounds(grid_x, grid_y, grid_z, radius):
	var chunks_to_update = {}
	var bounds_pad = radius + 1
	
	var min_chunk_x = floori(float(grid_x - bounds_pad) / chunk_size)
	var max_chunk_x = floori(float(grid_x + bounds_pad) / chunk_size)
	var min_chunk_y = floori(float(grid_y - bounds_pad) / chunk_size)
	var max_chunk_y = floori(float(grid_y + bounds_pad) / chunk_size)
	var min_chunk_z = floori(float(grid_z - bounds_pad) / chunk_size)
	var max_chunk_z = floori(float(grid_z + bounds_pad) / chunk_size)
	
	for cx in range(min_chunk_x, max_chunk_x + 1):
		for cy in range(min_chunk_y, max_chunk_y + 1):
			for cz in range(min_chunk_z, max_chunk_z + 1):
				if cx >= 0 and cx < chunk_count_x and \
				   cy >= 0 and cy < chunk_count_y and \
				   cz >= 0 and cz < chunk_count_z:
					var index = to_index(cx, cy, cz)
					if index >= 0 and index < chunks.size():
						var chunk = chunks[index]
						if is_instance_valid(chunk):
							chunks_to_update[index] = chunk

	for chunk in chunks_to_update.values():
		chunk.generate_mesh_from_density_grid()
