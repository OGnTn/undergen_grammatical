@tool
extends Control
## Visual Pipeline Graph Editor for UnderGenPipeline resources.
## Each node in the GraphEdit represents one UnderGenNode resource.
## Connections are serialized directly into the UnderGenPipeline resource.

var current_pipeline:
	set(value):
		current_pipeline = value
		_load_pipeline()

# --- Node Type Catalogue ---
# Maps display name -> { "class": ClassName, "color": Color, "description": String }
const NODE_CATALOGUE = {
	# ── Logic ────────────────────────────────────────────────────────────────
	"Grammar Expander": {
		"class": "UnderGenGrammarNode",
		"color": Color(0.3, 0.6, 0.9),
		"description": "Runs graph-grammar expansion on a LevelGrammarResource.\nConnect Logical Graph output → BSP Room Placer port 1.",
		"inputs":  [{"name": "Seed", "type": 0}],
		"outputs": [{"name": "Logical Graph", "type": 1}]
	},
	# ── Layout ────────────────────────────────────────────────────────────────
	"BSP Room Placer": {
		"class": "UnderGenBSPPlacerNode",
		"color": Color(0.2, 0.5, 0.8),
		"description": "Places rooms using Binary Space Partitioning.",
		"inputs": [{"name": "Seed", "type": 0}, {"name": "Logical Graph", "type": 1}],
		"outputs": [{"name": "Gen Context", "type": 2}]
	},
	"A* Carver": {
		"class": "UnderGenAStarCarverNode",
		"color": Color(0.6, 0.3, 0.8),
		"description": "Carves corridors using 3D A* pathfinding.",
		"inputs": [{"name": "Gen Context", "type": 2}],
		"outputs": [{"name": "Gen Context", "type": 2}]
	},
	"Bezier Carver": {
		"class": "UnderGenBezierCarverNode",
		"color": Color(0.8, 0.35, 0.6),
		"description": "Carves organic winding tunnels using bezier curves with wobble noise.",
		"inputs": [{"name": "Gen Context", "type": 2}],
		"outputs": [{"name": "Gen Context", "type": 2}]
	},
	"3D Noise Blend": {
		"class": "UnderGenNoiseNode",
		"color": Color(0.2, 0.6, 0.4),
		"description": "Blends 3D noise into the density grid.",
		"inputs": [{"name": "Gen Context", "type": 2}],
		"outputs": [{"name": "Gen Context", "type": 2}]
	},
	"Smooth Filter": {
		"class": "UnderGenSmoothNode",
		"color": Color(0.3, 0.5, 0.5),
		"description": "Applies a separable low-pass blur to the grid.",
		"inputs": [{"name": "Gen Context", "type": 2}],
		"outputs": [{"name": "Gen Context", "type": 2}]
	},
	"Vox Stamper": {
		"class": "UnderGenVoxStampNode",
		"color": Color(0.7, 0.5, 0.2),
		"description": "Stamps .vox files and extracts spawn markers.",
		"inputs": [{"name": "Gen Context", "type": 2}],
		"outputs": [{"name": "Gen Context", "type": 2}]
	},
	"Material Stamper": {
		"class": "UnderGenMaterialStamperNode",
		"color": Color(0.5, 0.65, 0.3),
		"description": "Stamps material IDs onto voxels based on zone-to-material mapping.",
		"inputs": [{"name": "Gen Context", "type": 2}],
		"outputs": [{"name": "Gen Context", "type": 2}]
	},
	"Surface Sampler": {
		"class": "UnderGenSurfaceSamplerNode",
		"color": Color(0.8, 0.6, 0.2),
		"description": "Samples surface voxels into a PointSet.",
		"inputs": [{"name": "Gen Context", "type": 2}],
		"outputs": [{"name": "Gen Context", "type": 2}, {"name": "Point Set", "type": 3}]
	},
	"Point Filter": {
		"class": "UnderGenPointFilterNode",
		"color": Color(0.5, 0.7, 0.3),
		"description": "Filters a PointSet by zone, slope, density, and spacing.",
		"inputs": [{"name": "Point Set", "type": 3}],
		"outputs": [{"name": "Point Set", "type": 3}]
	},
	"Marching Cubes Mesher": {
		"class": "UnderGenMesherNode",
		"color": Color(0.8, 0.2, 0.2),
		"description": "Terminal: Spawns MCChunk mesh nodes from the grid.",
		"inputs": [{"name": "Gen Context", "type": 2}],
		"outputs": []
	},
	"Scene Spawner": {
		"class": "UnderGenSceneSpawnerNode",
		"color": Color(0.8, 0.4, 0.1),
		"description": "Terminal: Spawns PackedScene instances at PointSet positions.",
		"inputs": [{"name": "Point Set", "type": 3}],
		"outputs": []
	},
}

# Port type colors (matches colour used in port rendering)
const PORT_COLORS = {
	0: Color(0.9, 0.9, 0.9), # Seed / Int
	1: Color(0.4, 0.8, 1.0), # Logical Graph (blue)
	2: Color(0.2, 0.9, 0.5), # Gen Context (green)
	3: Color(1.0, 0.85, 0.2), # PointSet (yellow)
}

# --- UI Components ---
var graph_edit: GraphEdit
var add_node_popup: PopupMenu
var _popup_open_position: Vector2 = Vector2.ZERO
var _node_counter: int = 0
var _grammar_pick_dlg: EditorFileDialog = null  # reused for grammar path picking
var _pending_grammar_node_res = null            # node waiting for grammar path

func _ready():
	_build_ui()

func _build_ui():
	set_anchors_preset(Control.PRESET_FULL_RECT)
	size_flags_horizontal = Control.SIZE_EXPAND_FILL
	size_flags_vertical = Control.SIZE_EXPAND_FILL

	var vbox = VBoxContainer.new()
	vbox.set_anchors_preset(Control.PRESET_FULL_RECT)
	vbox.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	vbox.size_flags_vertical = Control.SIZE_EXPAND_FILL
	add_child(vbox)

	# Toolbar
	var toolbar = HBoxContainer.new()
	vbox.add_child(toolbar)

	var new_btn = Button.new()
	new_btn.text = "✦ New Pipeline"
	new_btn.tooltip_text = "Create a new empty UnderGenPipeline resource."
	new_btn.pressed.connect(_on_new_pipeline)
	toolbar.add_child(new_btn)

	var load_btn = Button.new()
	load_btn.text = "📂 Load Pipeline"
	load_btn.tooltip_text = "Open an existing UnderGenPipeline resource from disk."
	load_btn.pressed.connect(_on_load_pipeline)
	toolbar.add_child(load_btn)

	var save_btn = Button.new()
	save_btn.text = "💾 Save Pipeline"
	save_btn.tooltip_text = "Save the current pipeline resource to disk."
	save_btn.pressed.connect(_on_save_pipeline)
	toolbar.add_child(save_btn)

	var separator = VSeparator.new()
	toolbar.add_child(separator)

	var hint = Label.new()
	hint.text = "Right-click on the canvas to add nodes. Click a node to inspect it."
	hint.add_theme_color_override("font_color", Color(0.7, 0.7, 0.7))
	toolbar.add_child(hint)

	# GraphEdit
	graph_edit = GraphEdit.new()
	graph_edit.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	graph_edit.size_flags_vertical = Control.SIZE_EXPAND_FILL
	graph_edit.right_disconnects = true
	graph_edit.connection_request.connect(_on_connection_request)
	graph_edit.disconnection_request.connect(_on_disconnection_request)
	graph_edit.popup_request.connect(_on_popup_request)
	graph_edit.delete_nodes_request.connect(_on_delete_nodes_request)
	graph_edit.node_selected.connect(_on_node_selected)
	vbox.add_child(graph_edit)

	# Add-node context menu
	add_node_popup = PopupMenu.new()
	add_node_popup.name = "AddNodePopup"
	var idx = 0
	for node_name in NODE_CATALOGUE:
		add_node_popup.add_item(node_name, idx)
		add_node_popup.set_item_tooltip(idx, NODE_CATALOGUE[node_name]["description"])
		idx += 1
	add_node_popup.id_pressed.connect(_on_add_node_selected)
	add_child(add_node_popup)

# ─────────────────────────────────────────────
#   Pipeline Loading / Saving
# ─────────────────────────────────────────────

func _load_pipeline():
	if not graph_edit: return
	graph_edit.clear_connections()

	# Remove all existing graph nodes
	for child in graph_edit.get_children():
		if child is GraphNode:
			child.queue_free()

	if not current_pipeline: return

	# Rebuild nodes
	var pipeline_nodes: Array = current_pipeline.get_nodes()
	for n in pipeline_nodes:
		if n == null: continue
		_create_graph_node_for(n)

	# Rebuild connections
	var pipeline_conns: Array = current_pipeline.get_connections()
	for conn in pipeline_conns:
		graph_edit.connect_node(
			conn.get("from_node", ""),
			conn.get("from_port", 0),
			conn.get("to_node", ""),
			conn.get("to_port", 0)
		)

func _save_pipeline_state():
	if not current_pipeline: return

	# Sync editor positions back into nodes
	for child in graph_edit.get_children():
		if child is GraphNode:
			var node_res = child.get_meta("node_resource", null)
			if node_res:
				node_res.set_editor_position(child.position_offset)

	# Sync connections
	current_pipeline.set_connections([]) # Clear
	var conns = graph_edit.get_connection_list()
	for conn in conns:
		current_pipeline.add_connection(conn["from_node"], conn["from_port"], conn["to_node"], conn["to_port"])

# ─────────────────────────────────────────────
#   Graph Node Creation (UI)
# ─────────────────────────────────────────────

func _create_graph_node_for(node_res) -> GraphNode:
	var class_name_str = node_res.get_class()

	# Find catalogue entry for display info
	var catalogue_key = ""
	for key in NODE_CATALOGUE:
		if NODE_CATALOGUE[key]["class"] == class_name_str:
			catalogue_key = key
			break

	var node_color = Color(0.3, 0.3, 0.3)
	var input_ports = []
	var output_ports = []

	if catalogue_key != "":
		node_color = NODE_CATALOGUE[catalogue_key]["color"]
		input_ports = NODE_CATALOGUE[catalogue_key].get("inputs", [])
		output_ports = NODE_CATALOGUE[catalogue_key].get("outputs", [])

	var gn = GraphNode.new()
	gn.title = catalogue_key if catalogue_key != "" else class_name_str
	gn.name = node_res.get_name() if node_res.get_name() != "" else "node_%d" % _node_counter
	_node_counter += 1
	gn.position_offset = node_res.get_editor_position()
	gn.set_meta("node_resource", node_res)

	# Color band
	gn.add_theme_color_override("title_color", node_color.lightened(0.4))

	# Build port rows
	var max_ports = max(input_ports.size(), output_ports.size())
	for i in range(max(1, max_ports)):
		var row = HBoxContainer.new()
		row.size_flags_horizontal = Control.SIZE_EXPAND_FILL

		var has_input = i < input_ports.size()
		var has_output = i < output_ports.size()

		# Input label
		var in_lbl = Label.new()
		if has_input:
			in_lbl.text = input_ports[i]["name"]
		row.add_child(in_lbl)

		# Spacer
		var sp = Control.new()
		sp.size_flags_horizontal = Control.SIZE_EXPAND_FILL
		row.add_child(sp)

		# Output label
		var out_lbl = Label.new()
		if has_output:
			out_lbl.text = output_ports[i]["name"]
			out_lbl.horizontal_alignment = HORIZONTAL_ALIGNMENT_RIGHT
		row.add_child(out_lbl)

		gn.add_child(row)

		# Register ports
		if has_input:
			var pt = input_ports[i]["type"]
			gn.set_slot(i, true, pt, PORT_COLORS.get(pt, Color.WHITE), has_output, output_ports[i]["type"] if has_output else 0, PORT_COLORS.get(output_ports[i]["type"] if has_output else 0, Color.WHITE))
		elif has_output:
			var pt = output_ports[i]["type"]
			gn.set_slot(i, false, 0, Color.WHITE, true, pt, PORT_COLORS.get(pt, Color.WHITE))

	graph_edit.add_child(gn)
	return gn

# ─────────────────────────────────────────────
#   Signals / Events
# ─────────────────────────────────────────────

func _on_connection_request(from_node: StringName, from_port: int, to_node: StringName, to_port: int):
	graph_edit.connect_node(from_node, from_port, to_node, to_port)
	if current_pipeline:
		current_pipeline.add_connection(from_node, from_port, to_node, to_port)

func _on_disconnection_request(from_node: StringName, from_port: int, to_node: StringName, to_port: int):
	graph_edit.disconnect_node(from_node, from_port, to_node, to_port)
	if current_pipeline:
		current_pipeline.remove_connection(from_node, from_port, to_node, to_port)

func _on_popup_request(at_position: Vector2):
	_popup_open_position = at_position + graph_edit.global_position
	add_node_popup.popup(Rect2(_popup_open_position, Vector2.ZERO))

func _on_add_node_selected(id: int):
	if not current_pipeline:
		printerr("PipelineEditor: No pipeline loaded. Create or assign one first.")
		return

	var node_names = NODE_CATALOGUE.keys()
	if id >= node_names.size(): return

	var key = node_names[id]
	var class_name_str: String = NODE_CATALOGUE[key]["class"]

	# Instantiate the C++ resource via ClassDB
	var node_res = ClassDB.instantiate(class_name_str)
	if node_res == null:
		printerr("PipelineEditor: Could not instantiate class '", class_name_str, "'. Check GDExtension is loaded.")
		return

	var graph_pos = ((_popup_open_position - graph_edit.global_position) + graph_edit.scroll_offset) / graph_edit.zoom
	node_res.set_editor_position(graph_pos)
	node_res.set_name(key.replace(" ", "_") + "_%d" % _node_counter)

	current_pipeline.add_pipeline_node(node_res)
	_create_graph_node_for(node_res)

	# Grammar Expander: immediately prompt for grammar resource path.
	if class_name_str == "UnderGenGrammarNode":
		_prompt_grammar_path_for(node_res)

func _on_delete_nodes_request(nodes: Array):
	for node_name in nodes:
		var gn = graph_edit.get_node_or_null(NodePath(node_name))
		if gn:
			var node_res = gn.get_meta("node_resource", null)
			if node_res and current_pipeline:
				current_pipeline.remove_pipeline_node(node_res)
			gn.queue_free()
	# Remove all connections to/from deleted nodes
	if current_pipeline:
		_save_pipeline_state()

func _on_node_selected(node: Node):
	# Forward selected resource to the Godot Inspector
	if node.has_meta("node_resource"):
		var res = node.get_meta("node_resource")
		EditorInterface.get_inspector().edit(res)

func _prompt_grammar_path_for(node_res) -> void:
	if not _grammar_pick_dlg:
		_grammar_pick_dlg = EditorFileDialog.new()
		_grammar_pick_dlg.access    = EditorFileDialog.ACCESS_RESOURCES
		_grammar_pick_dlg.file_mode = EditorFileDialog.FILE_MODE_OPEN_FILE
		_grammar_pick_dlg.add_filter("*.tres;*.res", "Grammar Resource")
		_grammar_pick_dlg.title     = "Select Grammar Resource for this node"
		add_child(_grammar_pick_dlg)
	for conn in _grammar_pick_dlg.file_selected.get_connections():
		_grammar_pick_dlg.file_selected.disconnect(conn.callable)
	_grammar_pick_dlg.file_selected.connect(func(path: String):
		node_res.set("grammar_resource_path", path)
		# Refresh graph node title to show grammar filename
		for child in graph_edit.get_children():
			if child is GraphNode and child.get_meta("node_resource", null) == node_res:
				child.title = "Grammar Expander\n" + path.get_file()
				break
		print("PipelineEditor: Grammar path set → ", path))
	_grammar_pick_dlg.popup_centered_ratio(0.65)

func _on_new_pipeline():
	var new_pipeline = ClassDB.instantiate("UnderGenPipeline")
	current_pipeline = new_pipeline

func _on_load_pipeline():
	var dlg = EditorFileDialog.new()
	dlg.access    = EditorFileDialog.ACCESS_RESOURCES
	dlg.file_mode = EditorFileDialog.FILE_MODE_OPEN_FILE
	dlg.add_filter("*.tres;*.res", "Godot Resource")
	dlg.title     = "Open Pipeline Resource"
	dlg.file_selected.connect(func(path: String):
		var res = load(path)
		if res != null and res.get_class() == "UnderGenPipeline":
			current_pipeline = res
		else:
			printerr("PipelineEditor: '", path, "' is not an UnderGenPipeline.")
		dlg.queue_free())
	dlg.canceled.connect(func(): dlg.queue_free())
	add_child(dlg)
	dlg.popup_centered_ratio(0.65)

func _on_save_pipeline():
	if not current_pipeline:
		printerr("PipelineEditor: Nothing to save.")
		return
	_save_pipeline_state()
	var path = current_pipeline.resource_path
	if path.is_empty():
		# Open a file dialog to choose save path
		var dialog = EditorFileDialog.new()
		dialog.access = EditorFileDialog.ACCESS_RESOURCES
		dialog.file_mode = EditorFileDialog.FILE_MODE_SAVE_FILE
		dialog.add_filter("*.tres", "Godot Resource")
		dialog.file_selected.connect(func(p): ResourceSaver.save(current_pipeline, p))
		add_child(dialog)
		dialog.popup_centered_ratio(0.6)
	else:
		ResourceSaver.save(current_pipeline, path)
		print("PipelineEditor: Saved to ", path)
