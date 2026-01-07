@tool
extends Node3D

var world: World
var multimeshes: Dictionary = {} # ZoneID (int) -> MultiMeshInstance3D

func _ready():
	# If added manually, try to find world
	if not world:
		world = get_parent() as World

func visualize():
	clear()
	
	if not world:
		print("ZoneVisualizer: No World assigned.")
		return
	if not world.density_grid_resource:
		print("ZoneVisualizer: No Density Grid.")
		return
		
	var grid = world.density_grid_resource
	
	print("ZoneVisualizer: Fetching optimized zone data from C++...")
	var zone_buckets = grid.get_surface_zones()
	
	print("ZoneVisualizer: Analyzed surface. Zones found: ", zone_buckets.size())

	# Create MultiMeshes
	for zone_id in zone_buckets:
		var points = zone_buckets[zone_id]
		if points.is_empty(): continue
		
		var mm_inst = MultiMeshInstance3D.new()
		var mm = MultiMesh.new()
		
		mm.transform_format = MultiMesh.TRANSFORM_3D
		mm.mesh = BoxMesh.new()
		mm.mesh.size = Vector3(world.cube_size * 0.5, world.cube_size * 0.5, world.cube_size * 0.5) # Half size boxes
		
		# Material with Color
		var mat = StandardMaterial3D.new()
		mat.albedo_color = _get_color_for_id(zone_id)
		mat.shading_mode = BaseMaterial3D.SHADING_MODE_UNSHADED # Bright flat color
		mm.mesh.surface_set_material(0, mat)
		
		mm.mesh.surface_set_material(0, mat)
		
		# Efficient C++ population
		grid.populate_multimesh(mm, points, world.cube_size)
			
		mm_inst.multimesh = mm
		mm_inst.name = "Zone_" + str(zone_id)
		add_child(mm_inst)
		multimeshes[zone_id] = mm_inst

		# Create Label at Centroid
		var raw_centroid = Vector3.ZERO
		for p in points: raw_centroid += p
		raw_centroid /= points.size()
		var centroid = (raw_centroid + Vector3(0.5, 0.5, 0.5)) * world.cube_size
		
		var label = Label3D.new()
		var zone_name = grid.get_zone_name_by_id(zone_id)
		if zone_name == "": zone_name = "Zone " + str(zone_id)
		
		label.text = zone_name
		label.billboard = BaseMaterial3D.BILLBOARD_ENABLED
		label.no_depth_test = true
		label.fixed_size = true
		label.font_size = 32
		label.outline_render_priority = 10
		label.modulate = Color.WHITE
		label.outline_modulate = Color.BLACK
		label.position = centroid + Vector3(0, 2.0, 0) # Slightly above
		add_child(label)
		
	print("ZoneVisualizer: visualization complete. Zones found: ", zone_buckets.size())

func clear():
	for child in get_children():
		child.queue_free()
	multimeshes.clear()

func _get_color_for_id(id: int) -> Color:
	if id == 0: return Color.GRAY # Default
	
	# Procedural random color based on ID seed
	var rng = RandomNumberGenerator.new()
	rng.seed = id * 12345
	
	return Color.from_hsv(rng.randf(), 0.8, 1.0)
