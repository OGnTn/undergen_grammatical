extends Node3D

@export var world: World
@export var poisson_disk_sampler: PoissonDiskSampler
# @export var spawn_list: Array[WorldSpawnableObject] # MOVED TO LEVEL RESOURCE
@export var spawn_parent: Node3D

# Map scene paths to PackedScenes
var scene_cache: Dictionary[String, PackedScene] = {}
# Map scene paths back to their Resource (for data access during spawn)
var resource_cache: Dictionary[String, WorldSpawnableObject] = {}

signal finished_spawning

const TIME_BUDGET_MS = 4
var spawn_queue: Array[Dictionary] = []
var active_spawn_source: Array[WorldSpawnableObject] = []

@onready
var test_cylinder_scene = load("res://entities/test_cylinder/test_cylinder.tscn")

func _ready():
	set_process(false)

func prepare_spawn_list(p_spawn_list: Array[WorldSpawnableObject]):
	active_spawn_source = p_spawn_list
	print("DEBUG: prepare_spawn_list called. Processing %d items." % active_spawn_source.size())
	for i in range(active_spawn_source.size()):
		var spawn_res = active_spawn_source[i]
		
		if not spawn_res or not spawn_res.scene: 
			print("DEBUG: Item %d is null or missing scene, skipping" % i)
			continue
		
		# Use the scene's resource path as the unique key
		var path_key = spawn_res.scene.resource_path
		
		# 1. Cache the scene (if spawn_res.scene is already a PackedScene, we just store it)
		if not scene_cache.has(path_key):
			scene_cache[path_key] = spawn_res.scene
		
		resource_cache[path_key] = spawn_res
		
		# 2. Add to ObjectSpawner whitelist
		if %ObjectSpawner:
			%ObjectSpawner.add_spawnable_scene(path_key)
		else:
			printerr("ERROR: %ObjectSpawner is missing!")
			
	%ObjectSpawner.spawn_function = spawn_at_location
	
	if multiplayer.is_server():
		if has_node("/root/Network"):
			Network.server_tell_im_ready_to_spawn.rpc()
			
	print("DEBUG: prepare_spawn_list COMPLETE")

func populate_level():
	if !multiplayer.is_server():
		return

	print("DEBUG: Server detected. Generating spawn points...")
	
	if resource_cache.is_empty():
		printerr("Populator: populate_level called but resource_cache is empty! Make sure World called prepare_spawn_list first.")
		return

	# Generate candidate points on the surface
	spawn_queue.clear()
	
	# --- 1. Vox Fixed Spawns ---
	var vox_spawns = world.density_grid_resource.get_vox_spawns()
	var spawn_map = world.level_resource.vox_spawn_map
	
	if not vox_spawns.is_empty():
		print("DEBUG: Found %d fixed Vox spawns" % vox_spawns.size())
		for spawn in vox_spawns:
			var idx = spawn["palette_index"]
			var grid_pos = spawn["position"]
			
			if spawn_map.has(idx):
				var spawn_res = spawn_map[idx]
				if spawn_res and spawn_res.scene:
					var path_key = spawn_res.scene.resource_path
					
					# Ensure it's in our local cache so spawining doesn't fail
					if not scene_cache.has(path_key):
						scene_cache[path_key] = spawn_res.scene
						resource_cache[path_key] = spawn_res
						if %ObjectSpawner: %ObjectSpawner.add_spawnable_scene(path_key)
					
					spawn_queue.append({
						"scene_path": path_key,
						"point": grid_pos, # Used by spawn_at_location to set position * cube_size
						# Fixed spawns usually don't need collision/rarity checks
					})

	# --- 2. Procedural Surface Spawns ---
	var points = poisson_disk_sampler.poisson_disk_sampling(world.surface_normals.keys(), 1)
 

	for point in points:
		var point_normal = world.surface_normals[point]
		
		if point_normal == Vector3.ZERO: continue
		
		var true_pos = world.get_true_surface_position(point)
		#var t = test_cylinder_scene.instantiate()
		#t.global_position = true_pos
		#add_child(t)
		# --- NEW: ZONE QUERY LOGIC ---
		# 1. Convert World Position to Grid Coordinate
		var lookup_pos = true_pos + (point_normal * (world.cube_size * 0.5))
		var grid_pos = Vector3i(lookup_pos / world.cube_size)
		#var grid_pos = Vector3i(true_pos / world.cube_size)
		
		# 2. Ask C++ for the Zone ID and convert to Name
		# Note: get_zone_at and get_zone_name_by_id must be bound in your C++ DensityGrid
		var zone_id_int = world.density_grid_resource.get_zone_at(grid_pos)
		var zone_name = world.density_grid_resource.get_zone_name_by_id(zone_id_int)
		# 3. Filter candidates based on Slope AND Zone
		var valid_candidates = active_spawn_source.filter(func(s): 
			# A. Check Slope (Existing)
			var slope_ok = point_normal.y >= s.normal_y_min and point_normal.y <= s.normal_y_max
			if not slope_ok: 
				print("Slope not okay")
				return false
			
			# B. Check Zone Name (New)
			# If required_zone is set, it MUST match the current zone_name
			if s.required_zones != [""] and zone_name not in s.required_zones:
				print("Not required zone")
				return false
				
			return true
		)
		
		if valid_candidates.is_empty(): 
			continue
		
		# 4. Shuffle candidates so high-index items don't get priority
		valid_candidates.shuffle()
		
		# 5. Iterate candidates and roll for rarity/collision
		for candidate in valid_candidates:
			# Rarity Check
			if randf() > candidate.spawn_chance:
				continue 
			
			# Pre-calculate transform
			var target_transform = calculate_transform(true_pos, point_normal, candidate)
			
			# Collision Check
			var is_valid_spot = true
			if candidate.check_collision:
				is_valid_spot = is_space_clear_voxel(target_transform, candidate)
			
			if is_valid_spot:
				spawn_queue.append({
					"scene_path": candidate.scene.resource_path, # Use consistent key
					"point": true_pos,
					"point_normal": point_normal,
					"final_transform": target_transform
				})
				break # Found an object for this point, move to next point

	if !spawn_queue.is_empty():
		print("DEBUG: Queued %d objects. Enabling process..." % spawn_queue.size())
		set_process(true)
	else:
		_on_spawning_finished()

func _process(_delta):
	if spawn_queue.is_empty():
		_on_spawning_finished()
		return

	var start_time = Time.get_ticks_msec()

	while !spawn_queue.is_empty():
		var data = spawn_queue.front()
		var path_key = data["scene_path"]
		
		if not resource_cache.has(path_key):
			spawn_queue.pop_front() 
			continue

		# Delegate actual instantiation to the Spawner
		%ObjectSpawner.spawn(data)
		
		spawn_queue.pop_front()

		if Time.get_ticks_msec() - start_time > TIME_BUDGET_MS:
			break
			
	if spawn_queue.is_empty():
		_on_spawning_finished()

func calculate_transform(pos: Vector3, normal: Vector3, settings: WorldSpawnableObject) -> Transform3D:
	var t = Transform3D()
	
	# 1. Position
	var world_pos = pos * world.cube_size 
	t.origin = world_pos 
	
	# 2. Alignment
	if settings.align_to_surface_normal:
		t = align_with_y(t, normal)
	
	# 3. Apply Local Offset (Rotated by alignment)
	t.origin += t.basis * settings.position_offset

	# 4. Random Yaw (Applied locally to preserve normal alignment)
	if settings.random_y_rotation:
		t = t.rotated_local(Vector3.UP, randf() * TAU)
		
	return t

func is_space_clear_voxel(target_transform: Transform3D, settings: WorldSpawnableObject) -> bool:
	print("Checking if space is clear")
	var density_grid = world.density_grid_resource
	var threshold = world.surface_threshold
	var cube_size = world.cube_size
	
	var bounds = settings.scan_bounds
	if !settings.align_to_surface_normal: # Assuming gravity object
		bounds.end.y -= 0.2 
	
	var start = bounds.position
	var end = bounds.end
	var step = cube_size
	
	var epsilon = 0.01 
	
	var x = start.x
	while x <= end.x + epsilon:
		var y = start.y
		while y <= end.y + epsilon:
			var z = start.z
			while z <= end.z + epsilon:
				
				# Clamp the check to the exact bounds so we don't overshoot
				var check_pos_local = Vector3(
					min(x, end.x), 
					min(y, end.y), 
					min(z, end.z)
				)
				
				var global_check_pos = target_transform * check_pos_local
				var grid_pos = (global_check_pos / cube_size).floor()
				
				# Safety check for grid bounds if needed, or rely on C++ clamping
				var density = density_grid.get_cell(Vector3i(grid_pos), 1.0)
				
				# Check for collision (Density > threshold means solid ground)
				# We want empty space (Density < threshold)
				if density > (threshold - 0.1):
					return false 

				z += step
			y += step
		x += step
		
	return true

func spawn_at_location(data):
	var path_key = data["scene_path"]
	var scene = scene_cache[path_key]
	var instance = scene.instantiate()
	instance.name = str(randi())
	
	if data.has("final_transform"):
		instance.transform = data["final_transform"]
	else:
		instance.position = data["point"] * world.cube_size
	
	return instance

func align_with_y(xform: Transform3D, new_y: Vector3) -> Transform3D:
	xform.basis.y = new_y
	
	# Fix for Z-facing walls:
	# If the normal (new_y) is parallel to the default Z axis (0,0,1),
	# the cross product will be zero. We must use a different axis (X) as a fallback.
	if abs(new_y.dot(Vector3(0, 0, 1))) > 0.99:
		# Parallel to Z: Use X-axis (1,0,0) to compute the Z-axis first
		xform.basis.z = Vector3(1, 0, 0).cross(new_y).normalized()
		xform.basis.x = xform.basis.y.cross(xform.basis.z).normalized()
	else:
		# Standard case: Use Z-axis (0,0,1) to compute the X-axis
		xform.basis.x = xform.basis.y.cross(Vector3(0, 0, 1)).normalized()
		xform.basis.z = xform.basis.x.cross(xform.basis.y).normalized()
		
	return xform

func _on_spawning_finished():
	if is_processing():
		set_process(false)
	
	# RPC with call_local handles both remote notification and local signal emission
	peer_finished_spawning.rpc()

@rpc("any_peer", "call_local", "reliable")
func peer_finished_spawning():
	emit_signal("finished_spawning")
