@tool
extends Control

@onready var filename_label = $HBoxContainer/Label

@onready var graph_edit = $GraphEdit
@export
var current_resource: LevelGraphResource:
	set(value):
		current_resource = value
		if current_resource:
			print("Editor selected: " + current_resource.resource_path)
			# Update a label if you have one
			if filename_label:
				filename_label.text = current_resource.resource_path
				print("Resource_name: " + current_resource.resource_path)
			# Trigger the load function to visualize nodes
			load_graph_from_resource()

func _ready():
	# Setup Context Menu
	var popup = graph_edit.get_menu_hbox()
	var add_room_btn = Button.new()
	add_room_btn.text = "Add Room"
	add_room_btn.pressed.connect(_add_node.bind("room"))
	popup.add_child(add_room_btn)
	
	var add_path_btn = Button.new()
	add_path_btn.text = "Add Path"
	add_path_btn.pressed.connect(_add_node.bind("path"))
	popup.add_child(add_path_btn)

	# Connect Graph Signals
	graph_edit.connection_request.connect(_on_connection_request)
	graph_edit.disconnection_request.connect(_on_disconnection_request)

func _add_node(type: String):
	var node_scene
	if type == "room":
		node_scene = preload("res://systems/generation/architecture/editor/nodes/room_node.tscn").instantiate()
	else:
		node_scene = preload("res://systems/generation/architecture/editor/nodes/path_node.tscn").instantiate()
		
	graph_edit.add_child(node_scene)
	node_scene.position_offset = (graph_edit.scroll_offset + Vector2(100, 100))
	return node_scene

func _on_connection_request(from, from_port, to, to_port):
	graph_edit.connect_node(from, from_port, to, to_port)

func _on_disconnection_request(from, from_port, to, to_port):
	graph_edit.disconnect_node(from, from_port, to, to_port)

func save_to_resource():
	if !current_resource: return
	
	var data = {"nodes": [], "edges": []}
	var connection_list = graph_edit.get_connection_list()
	
	# 1. Parse Rooms (Nodes)
	for child in graph_edit.get_children():
		if child.has_method("get_data"):
			# Only save ROOM nodes as "nodes"
			# Path nodes are effectively edges in our logic
			if child.title == "Room":
				data["nodes"].append(child.get_data())

	# 2. Parse Paths (Edges)
	# Logic: Find [Room A] -> [Path Node] -> [Room B]
	for child in graph_edit.get_children():
		if(child.name != "_connection_layer"):
			if child.title == "Path / Corridor":
				var path_node_name = child.name
				
				# Find who connects INTO this path (Source Room)
				var source_room = _find_connection_source(path_node_name, connection_list)
				# Find who connects OUT of this path (Target Room)
				var target_room = _find_connection_target(path_node_name, connection_list)
				
				if source_room and target_room:
					var edge_data = {
						"from": source_room.get_data()["id"],
						"to": target_room.get_data()["id"],
						"type": child.get_data()["type"] # "spider_tunnel"
					}
					data["edges"].append(edge_data)

	current_resource.graph_data = data
	ResourceSaver.save(current_resource, current_resource.resource_path)
	print("Saved Level Graph!")

# Helpers to trace the connections through the Path Node
func _find_connection_source(node_name, list):
	for conn in list:
		if conn["to_node"] == node_name:
			return graph_edit.get_node(NodePath(conn["from_node"]))
	return null

func _find_connection_target(node_name, list):
	for conn in list:
		if conn["from_node"] == node_name:
			return graph_edit.get_node(NodePath(conn["to_node"]))
	return null

func load_graph_from_resource():
	if !current_resource: return
	
	# 1. Clear existing nodes
	graph_edit.clear_connections()
	for child in graph_edit.get_children():
		if child is GraphNode:
			child.queue_free()
			
	# 2. Re-create nodes from data
	var data = current_resource.graph_data
	if not data.has("nodes"): return
	
	# Create Rooms
	var id_to_node_map = {}
	
	for node_data in data["nodes"]:
		var node = _add_node("room") # Update _add_node to return the instance!
		# Restore values
		node.get_node("TypeEdit").text = node_data["type"]
		node.get_node("Dimensions/MinDimensions/WidthEdit").value = node_data["min_size"].x
		node.get_node("Dimensions/MinDimensions/HeightEdit").value = node_data["min_size"].y
		node.get_node("Dimensions/MinDimensions/DepthEdit").value = node_data["min_size"].z
		node.get_node("Dimensions/MaxDimensions/WidthEdit").value = node_data["max_size"].x
		node.get_node("Dimensions/MaxDimensions/HeightEdit").value = node_data["max_size"].y
		node.get_node("Dimensions/MaxDimensions/DepthEdit").value = node_data["max_size"].z
		
		node.position_offset = node_data.get("editor_offset", Vector2.ZERO)
		
		# Store ID for connection mapping
		node.node_id = node_data["id"] 
		id_to_node_map[node.node_id] = node.name

	# 3. Re-create Paths (Edges)
	if not data.has("edges"): return
	
	for edge in data["edges"]:
		# Create Path Node
		var path_node = _add_node("path")
		path_node.get_node("TypeEdit").text = edge["type"]
		
		# Ideally, you saved position for paths too. 
		# If not, place them between rooms:
		path_node.position_offset = (graph_edit.get_node(NodePath(id_to_node_map[edge["from"]])).position_offset + graph_edit.get_node(NodePath(id_to_node_map[edge["to"]])).position_offset) / 2
		
		# Re-connect Graph Lines
		var room_from = id_to_node_map[edge["from"]]
		var room_to = id_to_node_map[edge["to"]]
		
		# Connect Room A -> Path
		graph_edit.connect_node(room_from, 0, path_node.name, 0)
		# Connect Path -> Room B
		graph_edit.connect_node(path_node.name, 0, room_to, 0)
