extends Node3D 

var players_spawned = false 

func _ready() -> void:
	get_window().title = str(OS.get_process_id())
	%PlayerSpawner.spawn_function = spawn_player_scene
	Network.connect("players_loaded_level", _on_players_loaded_level)
	#%World.call_deferred("generate_world")
	$World.generate_world()

func spawn_player_scene(data): 
	$LoadingScreen.visible = false
	#var p_scene: Node3D = load("res://scenes/player/player_scene.tscn").instantiate()
	var p_scene: Node3D = load("res://entities/player/player_scene3.tscn").instantiate()
	#p_scene.name = str(data["id"]) + "|" + str(OS.get_process_id())
	p_scene.name = str(data["id"])
	p_scene.global_position = data["init_pos"]
	p_scene.set_multiplayer_authority(data["id"])
	#p_scene.setup_player()
	
	return p_scene

func _on_players_loaded_level(): 
	print(str(multiplayer.get_unique_id())+ ": " + "Players loaded level") 
	spawn_player_scenes()

func spawn_player_scenes(): 
	if(players_spawned): 
		return 
	print(str(multiplayer.get_unique_id())+ ": " + "Spawning player scenes") 
	for p in Network.peers.keys(): 
		var data = {"id": p}
		#var player = %MultiplayerSpawner.spawn(data)
		
		
		var spawn_pos = %Generator.spawn_position 
		$DebugLayer/StatsPanel.start_pos = spawn_pos
		$DebugLayer/StatsPanel.end_pos = %Generator.end_position 
		#$DebugLayer/StatsPanel.player = p_scene
		var spawn_pos_randomized = spawn_pos + Vector3(randf_range(-5, 5), randf_range(-5, 5), randf_range(-5, 5)) 
		
		# Use the World's get_closest_surface_point which now reads from C++ DensityGrid
		var surface_pos = $World.get_closest_surface_point(spawn_pos_randomized.x, spawn_pos_randomized.y, spawn_pos_randomized.z)
		var init_pos = Vector3(surface_pos.x * %World.cube_size, surface_pos.y * %World.cube_size + %World.cube_size * 2, surface_pos.z * %World.cube_size)
		print(str(multiplayer.get_unique_id())+ ": " + "Spawning player scene for peer: " + str(p) + " at " + str(Vector3(surface_pos.x * %World.cube_size, surface_pos.y * %World.cube_size + %World.cube_size * 2, surface_pos.z * %World.cube_size)))
		data["init_pos"] = init_pos
		var player = %PlayerSpawner.spawn(data)
		#p_scene.global_position = Vector3(surface_pos.x * %World.cube_size, surface_pos.y * %World.cube_size + %World.cube_size * 2, surface_pos.z * %World.cube_size)
		#player.set_multiplayer_authority(p)
		#p_scene.set_multiplayer_authority(p)
		#$Players.add_child(p_scene)
	players_spawned = true 
