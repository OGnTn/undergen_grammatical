@tool
extends Control

var current_resource: LevelGrammarResource:
	set(value):
		current_resource = value
		_load_rule_list()

var h_split: HSplitContainer
var rule_list: ItemList
var graph_edit: GraphEdit
var toolbar: HBoxContainer

var current_rule_idx: int = -1
var edge_metadata = {} # Stores "from_to" -> type mapping

func _ready():
	# Build UI Procedurally
	set_anchors_preset(Control.PRESET_FULL_RECT)
	size_flags_horizontal = Control.SIZE_EXPAND_FILL
	size_flags_vertical = Control.SIZE_EXPAND_FILL
	
	var vbox = VBoxContainer.new()
	vbox.set_anchors_preset(Control.PRESET_FULL_RECT)
	vbox.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	vbox.size_flags_vertical = Control.SIZE_EXPAND_FILL
	add_child(vbox)
	
	toolbar = HBoxContainer.new()
	var add_rule_btn = Button.new()
	add_rule_btn.text = "Add Rule"
	add_rule_btn.pressed.connect(_on_add_rule)
	toolbar.add_child(add_rule_btn)
	
	var save_btn = Button.new()
	save_btn.text = "Save Resource"
	save_btn.pressed.connect(_on_save)
	toolbar.add_child(save_btn)
	
	vbox.add_child(toolbar)
	
	h_split = HSplitContainer.new()
	h_split.size_flags_vertical = Control.SIZE_EXPAND_FILL
	h_split.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	vbox.add_child(h_split)
	
	rule_list = ItemList.new()
	rule_list.custom_minimum_size = Vector2(200, 0)
	rule_list.size_flags_vertical = Control.SIZE_EXPAND_FILL
	rule_list.item_selected.connect(_on_rule_selected)
	h_split.add_child(rule_list)
	
	var graph_container = VBoxContainer.new()
	graph_container.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	graph_container.size_flags_vertical = Control.SIZE_EXPAND_FILL
	h_split.add_child(graph_container)
	
	# Graph Toolbar
	var graph_toolbar = HBoxContainer.new()
	graph_container.add_child(graph_toolbar)
	
	var lbl = Label.new()
	lbl.text = "Target Symbol (LHS):"
	lbl.tooltip_text = "The symbol in the graph that this rule will replace."
	graph_toolbar.add_child(lbl)
	
	var lhs_edit = LineEdit.new()
	lhs_edit.name = "LHSEdit"
	lhs_edit.custom_minimum_size = Vector2(100, 0)
	lhs_edit.placeholder_text = "e.g. Room"
	lhs_edit.text_changed.connect(_on_lhs_changed)
	graph_toolbar.add_child(lhs_edit)
	
	var prob_lbl = Label.new()
	prob_lbl.text = "Probability:"
	prob_lbl.tooltip_text = "Chance of this rule being picked if multiple rules match the same symbol."
	graph_toolbar.add_child(prob_lbl)
	var prob_spin = SpinBox.new()
	prob_spin.name = "ProbSpin"
	prob_spin.step = 0.1
	prob_spin.value = 1.0
	prob_spin.value_changed.connect(_on_prob_changed)
	graph_toolbar.add_child(prob_spin)
	
	var cond_lbl = Label.new()
	cond_lbl.text = "Condition:"
	cond_lbl.tooltip_text = "Rule only applies if this condition is true (e.g. keys > 2)."
	graph_toolbar.add_child(cond_lbl)
	var cond_edit = LineEdit.new()
	cond_edit.name = "CondEdit"
	cond_edit.custom_minimum_size = Vector2(120, 0)
	cond_edit.tooltip_text = "Format: VAR OP VALUE (e.g. keys > 2)"
	cond_edit.placeholder_text = "keys > 2"
	cond_edit.text_changed.connect(_on_cond_changed)
	graph_toolbar.add_child(cond_edit)

	var act_lbl = Label.new()
	act_lbl.text = "Action Effects:"
	act_lbl.tooltip_text = "Variables to change when this rule is applied."
	graph_toolbar.add_child(act_lbl)
	
	var act_edit = LineEdit.new()
	act_edit.name = "ActionsEdit"
	act_edit.custom_minimum_size = Vector2(120, 0)
	act_edit.tooltip_text = "JSON Dictionary (e.g. {\"keys\": 1} to add 1 key)"
	act_edit.placeholder_text = "{\"var\": 1}"
	act_edit.text_changed.connect(_on_actions_changed)
	graph_toolbar.add_child(act_edit)
	
	var edge_lbl = Label.new()
	edge_lbl.text = "New Edge Zone:"
	edge_lbl.tooltip_text = "Zone Type assigned to newly created connections (e.g. corridor, bridge, secret)."
	graph_toolbar.add_child(edge_lbl)
	
	var edge_edit = LineEdit.new()
	edge_edit.name = "EdgeTypeEdit"
	edge_edit.custom_minimum_size = Vector2(100, 0)
	edge_edit.text = "corridor"
	edge_edit.placeholder_text = "corridor"
	graph_toolbar.add_child(edge_edit)

	graph_edit = GraphEdit.new()
	graph_edit.size_flags_vertical = Control.SIZE_EXPAND_FILL
	graph_edit.right_disconnects = true
	# graph_edit.connection_request.connect...
	graph_edit.connection_request.connect(_on_connection_request)
	graph_edit.disconnection_request.connect(_on_disconnection_request)
	graph_edit.popup_request.connect(_on_popup_request)
	graph_edit.delete_nodes_request.connect(_on_delete_nodes_request) # Handle deletion via GraphEdit
	
	# Context Menu
	var popup_menu = PopupMenu.new()
	popup_menu.name = "ContextMenu"
	popup_menu.add_item("Add Room Node", 0)
	popup_menu.id_pressed.connect(_on_context_menu_item_selected)
	add_child(popup_menu)
	
	graph_container.add_child(graph_edit)


func _load_rule_list():
	rule_list.clear()
	if !current_resource: return
	
	for i in range(current_resource.rules.size()):
		var rule = current_resource.rules[i]
		rule_list.add_item(rule.rule_name + " (" + rule.lhs_symbol + ")")

func _on_add_rule():
	if !current_resource:
		printerr("Grammar Editor: current_resource is NULL")
		return
	
	print("Adding new rule...")
	var new_rule = GraphRule.new()
	new_rule.rule_name = "Rule " + str(current_resource.rules.size())
	
	var new_list = current_resource.rules.duplicate()
	new_list.append(new_rule)
	current_resource.rules = new_list
	
	print("Rule Added. Count: ", current_resource.rules.size())
	_load_rule_list()
	
	if current_resource.rules.size() > 0:
		var new_idx = current_resource.rules.size() - 1
		rule_list.select(new_idx)
		_on_rule_selected(new_idx)

func _on_rule_selected(idx):
	current_rule_idx = idx
	_load_graph_from_rule(current_resource.rules[idx])

func _load_graph_from_rule(rule: GraphRule):
	graph_edit.clear_connections()
	for c in graph_edit.get_children():
		if c is GraphNode: c.queue_free()
		
	# Reset entry/exit if not present (backward compatibility handled by export default)

		
	# Update Toolbar
	var lhs = find_child("LHSEdit", true, false)
	if lhs: lhs.text = rule.lhs_symbol
	var prob = find_child("ProbSpin", true, false)
	if prob: prob.value = rule.probability
	var cond = find_child("CondEdit", true, false)
	if cond: cond.text = rule.condition
	var act = find_child("ActionsEdit", true, false)
	if act: act.text = JSON.stringify(rule.actions)
	
	# Load Nodes
	var id_map = {}
	for n_data in rule.rhs_nodes:
		var g_node = _create_graph_node_ui(n_data["symbol"])
		g_node.name = n_data["id"]
		g_node.position_offset = n_data.get("editor_pos", Vector2(100, 100))
		
		# Set properties in UI slots
		g_node.get_node("Content/SymbolEdit").text = n_data["symbol"]
		g_node.get_node("Content/VoxEdit").text = n_data.get("vox_path", "")
		
		# Set Sizes
		if n_data.has("min_size"):
			var min_s = n_data["min_size"]
			g_node.get_node("Content/SizeBox/MinSizeEdit").text = "%d,%d,%d" % [min_s.x, min_s.y, min_s.z]
		else:
			g_node.get_node("Content/SizeBox/MinSizeEdit").text = ""

		if n_data.has("max_size"):
			var max_s = n_data["max_size"]
			g_node.get_node("Content/SizeBox/MaxSizeEdit").text = "%d,%d,%d" % [max_s.x, max_s.y, max_s.z]
		else:
			g_node.get_node("Content/SizeBox/MaxSizeEdit").text = ""
		
		# Set Constraints
		var constraints = n_data.get("constraints", {})
		if constraints.has("relative_to"):
			g_node.get_node("Content/ConstraintBox/RelativeEdit").text = constraints["relative_to"]
			if constraints.has("offset"):
				var off = constraints["offset"]
				g_node.get_node("Content/ConstraintBox/RelativeOffEdit").text = "%.1f,%.1f,%.1f" % [off.x, off.y, off.z]
		
	# CONNECT TO
		if n_data.has("connect_to"):
			var ct = n_data["connect_to"]
			g_node.get_node("Content/ConnectBox/ConnectTargetEdit").text = ct.get("target", "")
		
		if constraints.has("fixed_pos"):
			g_node.get_node("Content/ConstraintBox/FixedCheck").button_pressed = true
			var v = constraints["fixed_pos"]
			g_node.get_node("Content/ConstraintBox/FixedEdit").text = "%.1f,%.1f,%.1f" % [v.x, v.y, v.z]
			
		# Set Roles
		if rule.entry_node_id == n_data["id"]:
			g_node.get_node("Content/RoleBox/EntryCheck").button_pressed = true
		if rule.exit_node_id == n_data["id"]:
			g_node.get_node("Content/RoleBox/ExitCheck").button_pressed = true

		
		graph_edit.add_child(g_node)
		id_map[n_data["id"]] = g_node.name

	# Load Edges
	edge_metadata.clear()
	for e_data in rule.rhs_edges:
		if id_map.has(e_data["from"]) and id_map.has(e_data["to"]):
			var from_node = id_map[e_data["from"]]
			var to_node = id_map[e_data["to"]]
			graph_edit.connect_node(from_node, 0, to_node, 0)
			
			var key = from_node + "_" + to_node
			edge_metadata[key] = e_data.get("type", "corridor")

func _on_popup_request(at_position):
	var p = find_child("ContextMenu", true, false)
	if p:
		p.position = Vector2i(get_screen_position() + at_position)
		p.popup()

func _on_context_menu_item_selected(id):
	if id == 0: # Add RHS Node
		_add_graph_node_ui()

func _on_delete_nodes_request(nodes):
	# GraphEdit tells us which nodes to delete
	for node_name in nodes:
		var node = graph_edit.get_node(NodePath(node_name))
		if node:
			# Remove connections first? GraphEdit usually handles it or we call disconnect.
			# Just freeing the node is usually enough for UI, we save later.
			node.queue_free()
			# Ideally we should also remove connections involving this node to be clean
			var conns = graph_edit.get_connection_list()
			for c in conns:
				if c.from_node == node_name or c.to_node == node_name:
					graph_edit.disconnect_node(c.from_node, c.from_port, c.to_node, c.to_port)
					
					var key = c.from_node + "_" + c.to_node
					if edge_metadata.has(key): edge_metadata.erase(key)


func _create_graph_node_ui(symbol_text = "Symbol"):
	var gn = GraphNode.new()
	gn.title = "Room Node"
	
	gn.set_slot(0, true, 0, Color.WHITE, true, 0, Color.WHITE)
	
	var vbox = VBoxContainer.new()
	vbox.name = "Content"
	gn.add_child(vbox)
	
	var lbl = Label.new()
	lbl.text = "Room Type / Symbol:"
	vbox.add_child(lbl)
	
	# ROLES UI
	var role_box = HBoxContainer.new()
	role_box.name = "RoleBox"
	vbox.add_child(role_box)
	
	var entry_chk = CheckBox.new()
	entry_chk.name = "EntryCheck"
	entry_chk.text = "Input"
	entry_chk.tooltip_text = "Incoming edges connect here"
	role_box.add_child(entry_chk)
	
	var exit_chk = CheckBox.new()
	exit_chk.name = "ExitCheck"
	exit_chk.text = "Output"
	exit_chk.tooltip_text = "Outgoing edges start here"
	role_box.add_child(exit_chk)

	
	var edit = LineEdit.new()
	edit.name = "SymbolEdit"
	edit.text = symbol_text
	edit.placeholder_text = "e.g. Hallway"
	edit.custom_minimum_size = Vector2(80, 0)
	vbox.add_child(edit)
	
	var vox = LineEdit.new()
	vox.name = "VoxEdit"
	vox.placeholder_text = "res://path/to/mesh.vox"
	vox.tooltip_text = "Optional: Path to .vox file"
	vbox.add_child(vox)
	
	# SIZE UI
	var sz_lbl = Label.new()
	sz_lbl.text = "Size Range (Min - Max):"
	sz_lbl.add_theme_font_size_override("font_size", 10)
	vbox.add_child(sz_lbl)

	var size_hbox = HBoxContainer.new()
	size_hbox.name = "SizeBox"
	vbox.add_child(size_hbox)
	
	var min_edit = LineEdit.new()
	min_edit.name = "MinSizeEdit"
	min_edit.placeholder_text = "5,5,5"
	min_edit.tooltip_text = "Min Size (x,y,z)"
	min_edit.custom_minimum_size = Vector2(50, 0)
	size_hbox.add_child(min_edit)

	var max_edit = LineEdit.new()
	max_edit.name = "MaxSizeEdit"
	max_edit.placeholder_text = "10,10,10"
	max_edit.tooltip_text = "Max Size (x,y,z)"
	max_edit.custom_minimum_size = Vector2(50, 0)
	size_hbox.add_child(max_edit)

	# CONSTRAINTS UI
	var cbox = VBoxContainer.new()
	cbox.name = "ConstraintBox"
	vbox.add_child(cbox)
	
	var sep = HSeparator.new()
	cbox.add_child(sep)
	
	var clbl = Label.new()
	clbl.text = "Spatial Constraints"
	clbl.add_theme_font_size_override("font_size", 10)
	cbox.add_child(clbl)
	
	# CONNECT TO UI
	var conn_box = VBoxContainer.new()
	conn_box.name = "ConnectBox"
	vbox.add_child(conn_box)
	
	var clbl2 = Label.new()
	clbl2.text = "Loop / Snap Connection"
	clbl2.add_theme_font_size_override("font_size", 10)
	conn_box.add_child(clbl2)
	
	var ct_edit = LineEdit.new()
	ct_edit.name = "ConnectTargetEdit"
	ct_edit.placeholder_text = "Snap to Symbol..."
	conn_box.add_child(ct_edit)
	
	var rel_edit = LineEdit.new()
	rel_edit.name = "RelativeEdit"
	rel_edit.placeholder_text = "Align with ID..."
	rel_edit.tooltip_text = "ID of another node in this rule"
	cbox.add_child(rel_edit)
	
	var rel_off = LineEdit.new()
	rel_off.name = "RelativeOffEdit"
	rel_off.placeholder_text = "Offset x,y,z"
	rel_off.tooltip_text = "Offset relative to the aligned node"
	cbox.add_child(rel_off)

	
	var fix_chk = CheckBox.new()
	fix_chk.name = "FixedCheck"
	fix_chk.text = "Lock Position"
	fix_chk.tooltip_text = "Force this node to a specific world coordinate."
	cbox.add_child(fix_chk)
	
	var fix_edit = LineEdit.new()
	fix_edit.name = "FixedEdit"
	fix_edit.placeholder_text = "x,y,z"
	cbox.add_child(fix_edit)
	
	return gn


func _add_graph_node_ui():
	var gn = _create_graph_node_ui()
	graph_edit.add_child(gn)
	gn.position_offset = graph_edit.scroll_offset + Vector2(200, 200)

func _on_connection_request(from, from_port, to, to_port):
	graph_edit.connect_node(from, from_port, to, to_port)
	# Store Type
	var type_edit = find_child("EdgeTypeEdit", true, false)
	var type = "corridor"
	if type_edit: type = type_edit.text
	if type == "": type = "corridor"
	
	var key = from + "_" + to
	edge_metadata[key] = type
	print("Connected with Edge Type: ", type)

func _on_disconnection_request(from, from_port, to, to_port):
	graph_edit.disconnect_node(from, from_port, to, to_port)
	var key = from + "_" + to
	if edge_metadata.has(key): edge_metadata.erase(key)

func _on_lhs_changed(txt):
	if current_rule_idx < 0: return
	current_resource.rules[current_rule_idx].lhs_symbol = txt
	rule_list.set_item_text(current_rule_idx, current_resource.rules[current_rule_idx].rule_name + " (" + txt + ")")

func _on_prob_changed(val):
	if current_rule_idx < 0: return
	current_resource.rules[current_rule_idx].probability = val

func _on_cond_changed(txt):
	if current_rule_idx < 0: return
	current_resource.rules[current_rule_idx].condition = txt

func _on_actions_changed(txt):
	if current_rule_idx < 0: return
	# Simple JSON Check
	var p = JSON.parse_string(txt)
	if p == null:
		# Invalid JSON, maybe just ignore or change color?
		# find_child("ActionsEdit").modulate = Color(1,0,0)
		return
	if p is Dictionary:
		current_resource.rules[current_rule_idx].actions = p
		# find_child("ActionsEdit").modulate = Color(1,1,1)

func _on_save():
	if !current_resource or current_rule_idx < 0:
		if current_resource: ResourceSaver.save(current_resource)
		return
	
	# Serialize Current Graph to Rule
	var rule = current_resource.rules[current_rule_idx]
	rule.rhs_nodes.clear()
	rule.rhs_edges.clear()
	rule.entry_node_id = ""
	rule.exit_node_id = ""

	
	var connections = graph_edit.get_connection_list()
	
	for child in graph_edit.get_children():
		if child is GraphNode:
			var node_data = {
				"id": child.name,
				"symbol": child.get_node("Content/SymbolEdit").text,
				"vox_path": child.get_node("Content/VoxEdit").text,
				"editor_pos": child.position_offset,
				"constraints": {}
			}
			
			# SAVE SIZES
			var min_str = child.get_node("Content/SizeBox/MinSizeEdit").text
			var max_str = child.get_node("Content/SizeBox/MaxSizeEdit").text
			
			var min_parts = min_str.split(",")
			if min_parts.size() == 3:
				node_data["min_size"] = Vector3i(min_parts[0].to_int(), min_parts[1].to_int(), min_parts[2].to_int())
				
			var max_parts = max_str.split(",")
			if max_parts.size() == 3:
				node_data["max_size"] = Vector3i(max_parts[0].to_int(), max_parts[1].to_int(), max_parts[2].to_int())
			
			# SAVE CONSTRAINTS
			var rel = child.get_node("Content/ConstraintBox/RelativeEdit").text
			if rel != "":
				node_data["constraints"]["relative_to"] = rel
				var off_str = child.get_node("Content/ConstraintBox/RelativeOffEdit").text
				var o_parts = off_str.split(",")
				if o_parts.size() == 3:
					node_data["constraints"]["offset"] = Vector3(o_parts[0].to_float(), o_parts[1].to_float(), o_parts[2].to_float())
			
			# SAVE CONNECT TO
			var ct_target = child.get_node("Content/ConnectBox/ConnectTargetEdit").text
			if ct_target != "":
				node_data["connect_to"] = {
					"target": ct_target
				}
			
			if child.get_node("Content/ConstraintBox/FixedCheck").button_pressed:
				var vec_str = child.get_node("Content/ConstraintBox/FixedEdit").text
				var parts = vec_str.split(",")
				if parts.size() == 3:
					node_data["constraints"]["fixed_pos"] = Vector3(parts[0].to_float(), parts[1].to_float(), parts[2].to_float())
			
			# SAVE ROLES
			if child.get_node("Content/RoleBox/EntryCheck").button_pressed:
				rule.entry_node_id = node_data["id"]
			if child.get_node("Content/RoleBox/ExitCheck").button_pressed:
				rule.exit_node_id = node_data["id"]

			rule.rhs_nodes.append(node_data)

			
	for conn in connections:
		var key = conn["from_node"] + "_" + conn["to_node"]
		var type = edge_metadata.get(key, "corridor")
		
		rule.rhs_edges.append({
			"from": conn["from_node"],
			"to": conn["to_node"],
			"type": type
		})
		
	ResourceSaver.save(current_resource)
	print("Saved Grammar Resource")
