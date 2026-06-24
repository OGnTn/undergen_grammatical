@tool
extends Control
## Grammar Editor — self-contained two-tab workspace.
##
## Top bar:  [ Grammar: path ]  [ New ]  [ Open ]  [ Save ]
## Tab 0:    Room Palette  — define room types (the vocabulary)
## Tab 1:    Layout Rules  — build expansion rules (the structure)
##
## No FileSystem dock interaction required. Everything is created and
## saved from within this editor.

# ─── active resource ────────────────────────────────────────────────────────

var current_resource: LevelGrammarResource:
	set(value):
		current_resource = value
		_on_resource_changed()

# ─── top-level UI refs ──────────────────────────────────────────────────────

var _path_label:    Label          # shows current resource path
var _tabs:          TabContainer   # shown only when resource is loaded
var _splash:        Control        # shown when no resource is loaded

# Palette tab
var _pal_list:      ItemList
var _pal_insp:      ScrollContainer
var _pal_insp_box:  VBoxContainer
var _pal_sel:       int = -1
var _pal_file_dlg:  EditorFileDialog = null

# Rules tab
var _rules_list:    ItemList
var _rules_sel:     int = -1
var _edge_meta:     Dictionary = {}

var _state_flow:    HFlowContainer

var _lhs_opt:       OptionButton
var _prob_spin:     SpinBox
var _cond_chk:      CheckBox
var _cond_row:      HBoxContainer
var _cond_var_opt:  OptionButton
var _cond_op_opt:   OptionButton
var _cond_val_sp:   SpinBox
var _zone_opt:      OptionButton
var _actions_vbox:  VBoxContainer

var _graph:         GraphEdit
var _ctx_menu:      PopupMenu
var _ctx_pos:       Vector2
var _node_counter:  int = 0

# File dialogs (created on demand)
var _new_dlg:  EditorFileDialog = null
var _open_dlg: EditorFileDialog = null

const ZONE_LABELS := ["corridor", "bridge", "secret", "shortcut", "tunnel"]
const OP_LABELS   := ["<", ">", "<=", ">=", "==", "!="]

# ─────────────────────────────────────────────────────────────────────────────
#  BUILD UI
# ─────────────────────────────────────────────────────────────────────────────

func _ready():
	set_anchors_preset(Control.PRESET_FULL_RECT)
	size_flags_horizontal = Control.SIZE_EXPAND_FILL
	size_flags_vertical   = Control.SIZE_EXPAND_FILL

	var vbox = VBoxContainer.new()
	vbox.set_anchors_preset(Control.PRESET_FULL_RECT)
	vbox.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	vbox.size_flags_vertical   = Control.SIZE_EXPAND_FILL
	add_child(vbox)

	# ── Top resource bar ─────────────────────────────────────────────────────
	var top_bar = HBoxContainer.new()
	top_bar.add_theme_constant_override("separation", 6)
	vbox.add_child(top_bar)

	var grammar_lbl = Label.new()
	grammar_lbl.text = "Grammar:"
	top_bar.add_child(grammar_lbl)

	_path_label = Label.new()
	_path_label.text = "(no grammar loaded)"
	_path_label.modulate = Color(0.6, 0.6, 0.6)
	_path_label.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	_path_label.clip_text = true
	top_bar.add_child(_path_label)

	var new_btn = Button.new()
	new_btn.text = "✦ New"
	new_btn.tooltip_text = "Create a new grammar resource (.tres) inside your project."
	new_btn.pressed.connect(_on_new_grammar)
	top_bar.add_child(new_btn)

	var open_btn = Button.new()
	open_btn.text = "📂 Open"
	open_btn.tooltip_text = "Open an existing LevelGrammarResource (.tres) file."
	open_btn.pressed.connect(_on_open_grammar)
	top_bar.add_child(open_btn)

	var save_btn = Button.new()
	save_btn.text = "💾 Save"
	save_btn.tooltip_text = "Save the current grammar to disk."
	save_btn.pressed.connect(_on_save)
	top_bar.add_child(save_btn)

	top_bar.add_child(VSeparator.new())

	var import_json_btn = Button.new()
	import_json_btn.text = "⤓ Import JSON"
	import_json_btn.tooltip_text = "Load a .grammar.json file and convert to a grammar resource."
	import_json_btn.pressed.connect(_on_import_json)
	top_bar.add_child(import_json_btn)

	var export_json_btn = Button.new()
	export_json_btn.text = "⤒ Export JSON"
	export_json_btn.tooltip_text = "Export the current grammar as a .grammar.json file."
	export_json_btn.pressed.connect(_on_export_json)
	top_bar.add_child(export_json_btn)

	vbox.add_child(HSeparator.new())

	# ── Splash (no resource loaded) ──────────────────────────────────────────
	_splash = _build_splash()
	vbox.add_child(_splash)

	# ── Tab container (resource loaded) ──────────────────────────────────────
	_tabs = TabContainer.new()
	_tabs.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	_tabs.size_flags_vertical   = Control.SIZE_EXPAND_FILL
	_tabs.hide()
	vbox.add_child(_tabs)

	_tabs.add_child(_build_palette_tab())
	_tabs.add_child(_build_rules_tab())
	_tabs.set_tab_title(0, "🎨 Room Palette")
	_tabs.set_tab_title(1, "📐 Layout Rules")


func _build_splash() -> Control:
	var center = CenterContainer.new()
	center.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	center.size_flags_vertical   = Control.SIZE_EXPAND_FILL

	var vbox = VBoxContainer.new()
	vbox.alignment = BoxContainer.ALIGNMENT_CENTER
	center.add_child(vbox)

	var icon_lbl = Label.new()
	icon_lbl.text = "🗺"
	icon_lbl.add_theme_font_size_override("font_size", 48)
	icon_lbl.horizontal_alignment = HORIZONTAL_ALIGNMENT_CENTER
	vbox.add_child(icon_lbl)

	var title = Label.new()
	title.text = "No Grammar Loaded"
	title.add_theme_font_size_override("font_size", 20)
	title.horizontal_alignment = HORIZONTAL_ALIGNMENT_CENTER
	vbox.add_child(title)

	var sub = Label.new()
	sub.text = "Create a new grammar or open an existing one to get started."
	sub.modulate = Color(0.7, 0.7, 0.7)
	sub.horizontal_alignment = HORIZONTAL_ALIGNMENT_CENTER
	vbox.add_child(sub)

	vbox.add_child(HSeparator.new())

	var btn_row = HBoxContainer.new()
	btn_row.alignment = BoxContainer.ALIGNMENT_CENTER
	vbox.add_child(btn_row)

	var new_big = Button.new()
	new_big.text = "✦  New Grammar"
	new_big.custom_minimum_size = Vector2(160, 40)
	new_big.pressed.connect(_on_new_grammar)
	btn_row.add_child(new_big)

	var open_big = Button.new()
	open_big.text = "📂  Open Grammar"
	open_big.custom_minimum_size = Vector2(160, 40)
	open_big.pressed.connect(_on_open_grammar)
	btn_row.add_child(open_big)

	return center

# ─────────────────────────────────────────────────────────────────────────────
#  RESOURCE MANAGEMENT
# ─────────────────────────────────────────────────────────────────────────────

func _on_resource_changed():
	if current_resource:
		_path_label.text = current_resource.resource_path \
			if current_resource.resource_path != "" else "(unsaved)"
		_path_label.modulate = Color(0.9, 0.9, 0.9)
		_splash.hide()
		_tabs.show()
		_refresh_all()
	else:
		_path_label.text = "(no grammar loaded)"
		_path_label.modulate = Color(0.6, 0.6, 0.6)
		_splash.show()
		_tabs.hide()


func _on_new_grammar():
	if not _new_dlg:
		_new_dlg = EditorFileDialog.new()
		_new_dlg.access    = EditorFileDialog.ACCESS_RESOURCES
		_new_dlg.file_mode = EditorFileDialog.FILE_MODE_SAVE_FILE
		_new_dlg.add_filter("*.tres", "Godot Resource")
		add_child(_new_dlg)
	# Always reconnect — other paths (save-as) may have overridden the callback
	for conn in _new_dlg.file_selected.get_connections():
		_new_dlg.file_selected.disconnect(conn.callable)
	_new_dlg.title = "Create New Grammar Resource"
	_new_dlg.file_selected.connect(_create_grammar_at_path)
	_new_dlg.popup_centered_ratio(0.65)


func _create_grammar_at_path(path: String):
	var res = LevelGrammarResource.new()
	# Add a default "Start" room type so the axiom makes sense immediately.
	var start_rt = RoomType.new()
	start_rt.symbol = "Start"
	start_rt.color  = Color(0.3, 0.7, 0.5)
	res.room_types.append(start_rt)
	var err = ResourceSaver.save(res, path)
	if err != OK:
		printerr("GrammarEditor: Failed to save new grammar to ", path)
		return
	# Reload so resource_path is populated.
	current_resource = load(path) as LevelGrammarResource


func _on_open_grammar():
	if not _open_dlg:
		_open_dlg = EditorFileDialog.new()
		_open_dlg.access    = EditorFileDialog.ACCESS_RESOURCES
		_open_dlg.file_mode = EditorFileDialog.FILE_MODE_OPEN_FILE
		_open_dlg.add_filter("*.tres;*.res", "Godot Resource")
		_open_dlg.title     = "Open Grammar Resource"
		_open_dlg.file_selected.connect(_load_grammar_from_path)
		add_child(_open_dlg)
	_open_dlg.popup_centered_ratio(0.65)


func _load_grammar_from_path(path: String):
	var res = load(path)
	if res is LevelGrammarResource:
		current_resource = res
	else:
		printerr("GrammarEditor: '", path, "' is not a LevelGrammarResource.")


func _on_save():
	_save_current_rule()
	if not current_resource:
		printerr("GrammarEditor: No resource loaded.")
		return
	if current_resource.resource_path == "":
		# Never been saved — open save-as dialog for the current resource.
		if not _new_dlg:
			_new_dlg = EditorFileDialog.new()
			_new_dlg.access    = EditorFileDialog.ACCESS_RESOURCES
			_new_dlg.file_mode = EditorFileDialog.FILE_MODE_SAVE_FILE
			_new_dlg.add_filter("*.tres", "Godot Resource")
			add_child(_new_dlg)
		for conn in _new_dlg.file_selected.get_connections():
			_new_dlg.file_selected.disconnect(conn.callable)
		_new_dlg.title = "Save Grammar As…"
		_new_dlg.file_selected.connect(_save_current_as)
		_new_dlg.popup_centered_ratio(0.65)
		return
	ResourceSaver.save(current_resource)
	print("GrammarEditor: Saved → ", current_resource.resource_path)


func _save_current_as(path: String):
	if not current_resource: return
	var err = ResourceSaver.save(current_resource, path)
	if err != OK:
		printerr("GrammarEditor: Failed to save to ", path)
		return
	current_resource = load(path) as LevelGrammarResource
	_on_resource_changed()
	print("GrammarEditor: Saved → ", path)

# ── JSON Import / Export ─────────────────────────────────────────────────

var _json_import_dlg: EditorFileDialog = null
var _json_export_dlg: EditorFileDialog = null

func _on_import_json():
	if not _json_import_dlg:
		_json_import_dlg = EditorFileDialog.new()
		_json_import_dlg.access    = EditorFileDialog.ACCESS_RESOURCES
		_json_import_dlg.file_mode = EditorFileDialog.FILE_MODE_OPEN_FILE
		_json_import_dlg.add_filter("*.grammar.json;*.json", "Grammar JSON files")
		_json_import_dlg.title     = "Import Grammar from JSON"
		_json_import_dlg.file_selected.connect(_do_import_json)
		add_child(_json_import_dlg)
	_json_import_dlg.popup_centered_ratio(0.65)


func _do_import_json(path: String):
	# Use C++ LevelGrammarSpec for robust JSON parsing (avoids GDScript
	# typed-array serialization issues with ResourceSaver).
	var spec: RefCounted = ClassDB.instantiate("LevelGrammarSpec")
	if spec == null:
		printerr("GrammarEditor: LevelGrammarSpec not available. Is the GDExtension loaded?")
		return
	var err = spec.from_json(FileAccess.get_file_as_string(path))
	if err != OK:
		printerr("GrammarEditor: JSON parse error in ", path)
		return

	# Convert from C++ LevelGrammarSpec → GDScript LevelGrammarResource
	# (property-by-property copy avoids typed-array serialization bugs)
	var res = LevelGrammarResource.new()
	res.axiom = spec.axiom

	for sv in spec.state_variables:
		res.state_variables.append(sv)

	for rt in spec.room_types:
		var room = RoomType.new()
		room.symbol   = rt.symbol
		room.color    = rt.color
		room.weight   = rt.weight
		room.min_size = rt.min_size
		room.max_size = rt.max_size
		room.vox_path = rt.vox_path
		res.room_types.append(room)

	for rule in spec.rules:
		var r = GraphRule.new()
		r.rule_name      = rule.rule_name
		r.lhs_symbol     = rule.lhs_symbol
		r.probability    = rule.probability
		r.entry_node_id  = rule.entry_node_id
		r.exit_node_id   = rule.exit_node_id
		r.condition_var  = rule.condition_var
		r.condition_op   = rule.condition_op
		r.condition_val  = rule.condition_val
		# assign() is fine here — these won't be saved until the parent
		# resource is saved, and by then they're owned by GraphRule
		r.actions.assign(rule.actions)
		r.rhs_nodes.assign(rule.rhs_nodes)
		r.rhs_edges.assign(rule.rhs_edges)
		res.rules.append(r)

	current_resource = res
	print("GrammarEditor: Imported grammar from JSON → ", path)


func _on_export_json():
	if not current_resource:
		printerr("GrammarEditor: No grammar loaded to export.")
		return

	if not _json_export_dlg:
		_json_export_dlg = EditorFileDialog.new()
		_json_export_dlg.access    = EditorFileDialog.ACCESS_RESOURCES
		_json_export_dlg.file_mode = EditorFileDialog.FILE_MODE_SAVE_FILE
		_json_export_dlg.add_filter("*.grammar.json", "Grammar JSON files")
		_json_export_dlg.title     = "Export Grammar as JSON"
		_json_export_dlg.file_selected.connect(_do_export_json)
		add_child(_json_export_dlg)
	_json_export_dlg.popup_centered_ratio(0.65)


func _do_export_json(path: String):
	var d = _grammar_to_dictionary()
	var text = JSON.stringify(d, "\t")
	var file = FileAccess.open(path, FileAccess.WRITE)
	if file == null:
		printerr("GrammarEditor: Could not write to: ", path)
		return
	file.store_string(text)
	file.close()
	print("GrammarEditor: Exported grammar to JSON → ", path)


func _grammar_to_dictionary() -> Dictionary:
	var d: Dictionary = {}
	if current_resource == null:
		return d

	d["axiom"] = current_resource.axiom

	# State variables
	if not current_resource.state_variables.is_empty():
		d["state_variables"] = current_resource.state_variables.duplicate()

	# Room types
	var rt_arr: Array = []
	for rt: RoomType in current_resource.room_types:
		var r: Dictionary = {}
		r["symbol"]   = rt.symbol
		r["color"]    = [rt.color.r, rt.color.g, rt.color.b, rt.color.a]
		r["weight"]   = rt.weight
		r["min_size"] = [rt.min_size.x, rt.min_size.y, rt.min_size.z]
		r["max_size"] = [rt.max_size.x, rt.max_size.y, rt.max_size.z]
		if not rt.vox_path.is_empty():
			r["vox_path"] = rt.vox_path
		rt_arr.append(r)
	d["room_types"] = rt_arr

	# Rules
	var r_arr: Array = []
	for r: GraphRule in current_resource.rules:
		var rd: Dictionary = {}
		rd["rule_name"]    = r.rule_name
		rd["lhs_symbol"]   = r.lhs_symbol
		rd["probability"]  = r.probability
		if not r.entry_node_id.is_empty():
			rd["entry_node_id"] = r.entry_node_id
		if not r.exit_node_id.is_empty():
			rd["exit_node_id"]  = r.exit_node_id

		# Condition (only emit if var is set)
		if not r.condition_var.is_empty():
			rd["condition_var"] = r.condition_var
			rd["condition_op"]  = r.condition_op
			rd["condition_val"] = r.condition_val

		# Actions
		if not r.actions.is_empty():
			rd["actions"] = r.actions.duplicate()

		rd["rhs_nodes"] = r.rhs_nodes.duplicate(true)
		rd["rhs_edges"] = r.rhs_edges.duplicate(true)
		r_arr.append(rd)
	d["rules"] = r_arr

	return d

# ─────────────────────────────────────────────────────────────────────────────
#  PALETTE TAB
# ─────────────────────────────────────────────────────────────────────────────

func _build_palette_tab() -> Control:
	var root = HSplitContainer.new()
	root.name = "RoomPalette"
	root.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	root.size_flags_vertical   = Control.SIZE_EXPAND_FILL

	# Left: list
	var lbox = VBoxContainer.new()
	lbox.custom_minimum_size = Vector2(200, 0)
	root.add_child(lbox)

	var btn_bar = HBoxContainer.new()
	lbox.add_child(btn_bar)

	var add_btn = Button.new()
	add_btn.text = "+ Add Room Type"
	add_btn.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	add_btn.pressed.connect(_on_pal_add)
	btn_bar.add_child(add_btn)

	var del_btn = Button.new()
	del_btn.text = "🗑"
	del_btn.tooltip_text = "Delete selected room type"
	del_btn.pressed.connect(_on_pal_delete)
	btn_bar.add_child(del_btn)

	_pal_list = ItemList.new()
	_pal_list.size_flags_vertical = Control.SIZE_EXPAND_FILL
	_pal_list.item_selected.connect(_on_pal_selected)
	lbox.add_child(_pal_list)

	# Right: inspector
	_pal_insp = ScrollContainer.new()
	_pal_insp.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	_pal_insp.size_flags_vertical   = Control.SIZE_EXPAND_FILL
	root.add_child(_pal_insp)

	_pal_insp_box = VBoxContainer.new()
	_pal_insp_box.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	_pal_insp.add_child(_pal_insp_box)

	_pal_insp_box.add_child(_mk("Select a room type from the list to edit it."))
	return root


func _refresh_palette_list():
	_pal_list.clear()
	if not current_resource: return
	for rt: RoomType in current_resource.room_types:
		var idx = _pal_list.add_item(rt.symbol)
		_pal_list.set_item_icon(idx, _color_icon(rt.color))


func _on_pal_selected(idx: int):
	_pal_sel = idx
	if not current_resource or idx < 0 or idx >= current_resource.room_types.size(): return
	_load_pal_inspector(current_resource.room_types[idx])


func _on_pal_add():
	if not current_resource: return
	var rt = RoomType.new()
	rt.symbol = "Room%d" % current_resource.room_types.size()
	current_resource.room_types.append(rt)
	_refresh_palette_list()
	_refresh_lhs_options()
	_rebuild_ctx_menu()
	var new_idx = current_resource.room_types.size() - 1
	_pal_list.select(new_idx)
	_on_pal_selected(new_idx)


func _on_pal_delete():
	if not current_resource or _pal_sel < 0: return
	current_resource.room_types.remove_at(_pal_sel)
	_pal_sel = -1
	_refresh_palette_list()
	_refresh_lhs_options()
	_rebuild_ctx_menu()
	for c in _pal_insp_box.get_children(): c.queue_free()
	_pal_insp_box.add_child(_mk("Select a room type from the list to edit it."))


func _load_pal_inspector(rt: RoomType):
	for c in _pal_insp_box.get_children(): c.queue_free()

	var title = Label.new()
	title.text = "Room Type Properties"
	title.add_theme_font_size_override("font_size", 14)
	_pal_insp_box.add_child(title)
	_pal_insp_box.add_child(HSeparator.new())

	_pal_insp_box.add_child(_mk("Symbol (name used in rules):"))
	var sym_edit = LineEdit.new()
	sym_edit.text = rt.symbol
	sym_edit.text_changed.connect(func(t):
		rt.symbol = t
		_refresh_palette_list()
		_refresh_lhs_options()
		_rebuild_ctx_menu()
		_refresh_graph_node_colors())
	_pal_insp_box.add_child(sym_edit)

	_pal_insp_box.add_child(_mk("Color (shown in graph):"))
	var col_btn = ColorPickerButton.new()
	col_btn.color = rt.color
	col_btn.color_changed.connect(func(c):
		rt.color = c
		_refresh_palette_list()
		_refresh_graph_node_colors())
	_pal_insp_box.add_child(col_btn)

	_pal_insp_box.add_child(_mk("Weight (relative spawn frequency):"))
	var w_sp = SpinBox.new()
	w_sp.step = 0.1; w_sp.min_value = 0.0; w_sp.max_value = 100.0; w_sp.value = rt.weight
	w_sp.value_changed.connect(func(v): rt.weight = v)
	_pal_insp_box.add_child(w_sp)

	_pal_insp_box.add_child(HSeparator.new())
	_pal_insp_box.add_child(_mk("Min Size (x, y, z):"))
	_pal_insp_box.add_child(_vec3i_row(rt.min_size, func(v): rt.min_size = v))
	_pal_insp_box.add_child(_mk("Max Size (x, y, z):"))
	_pal_insp_box.add_child(_vec3i_row(rt.max_size, func(v): rt.max_size = v))

	_pal_insp_box.add_child(HSeparator.new())
	_pal_insp_box.add_child(_mk("Vox asset (optional):"))
	var vox_row = HBoxContainer.new()
	var vox_edit = LineEdit.new()
	vox_edit.text = rt.vox_path
	vox_edit.placeholder_text = "res://path/to/room.vox"
	vox_edit.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	vox_edit.text_changed.connect(func(t): rt.vox_path = t)
	vox_row.add_child(vox_edit)
	var browse = Button.new()
	browse.text = "Browse…"
	browse.pressed.connect(func(): _open_vox_picker(vox_edit, rt))
	vox_row.add_child(browse)
	_pal_insp_box.add_child(vox_row)


func _open_vox_picker(target_edit: LineEdit, rt: RoomType):
	if not _pal_file_dlg:
		_pal_file_dlg = EditorFileDialog.new()
		_pal_file_dlg.access    = EditorFileDialog.ACCESS_RESOURCES
		_pal_file_dlg.file_mode = EditorFileDialog.FILE_MODE_OPEN_FILE
		_pal_file_dlg.add_filter("*.vox", "MagicaVoxel files")
		add_child(_pal_file_dlg)
	for conn in _pal_file_dlg.file_selected.get_connections():
		_pal_file_dlg.file_selected.disconnect(conn.callable)
	_pal_file_dlg.file_selected.connect(func(path: String):
		target_edit.text = path; rt.vox_path = path)
	_pal_file_dlg.popup_centered_ratio(0.7)

# ─────────────────────────────────────────────────────────────────────────────
#  RULES TAB
# ─────────────────────────────────────────────────────────────────────────────

func _build_rules_tab() -> Control:
	var root = VBoxContainer.new()
	root.name = "LayoutRules"
	root.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	root.size_flags_vertical   = Control.SIZE_EXPAND_FILL

	# State variables bar
	var sv_row = HBoxContainer.new()
	sv_row.custom_minimum_size = Vector2(0, 32)
	root.add_child(sv_row)
	sv_row.add_child(_mk("State Vars:"))
	_state_flow = HFlowContainer.new()
	_state_flow.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	sv_row.add_child(_state_flow)
	var sv_add = Button.new(); sv_add.text = "+ Variable"
	sv_add.pressed.connect(_on_add_state_var)
	sv_row.add_child(sv_add)
	root.add_child(HSeparator.new())

	# Main split
	var split = HSplitContainer.new()
	split.size_flags_vertical   = Control.SIZE_EXPAND_FILL
	split.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	root.add_child(split)

	# Left: rule list
	var lbox = VBoxContainer.new()
	lbox.custom_minimum_size = Vector2(200, 0)
	split.add_child(lbox)
	var rbtn = HBoxContainer.new(); lbox.add_child(rbtn)
	var add_r = Button.new(); add_r.text = "+ Rule"
	add_r.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	add_r.pressed.connect(_on_add_rule); rbtn.add_child(add_r)
	var del_r = Button.new(); del_r.text = "🗑"
	del_r.tooltip_text = "Delete selected rule"
	del_r.pressed.connect(_on_delete_rule); rbtn.add_child(del_r)
	_rules_list = ItemList.new()
	_rules_list.size_flags_vertical = Control.SIZE_EXPAND_FILL
	_rules_list.item_selected.connect(_on_rule_selected)
	lbox.add_child(_rules_list)

	# Right: toolbar + graph
	var rbox = VBoxContainer.new()
	rbox.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	rbox.size_flags_vertical   = Control.SIZE_EXPAND_FILL
	split.add_child(rbox)
	rbox.add_child(_build_rule_toolbar())

	_graph = GraphEdit.new()
	_graph.size_flags_vertical = Control.SIZE_EXPAND_FILL
	_graph.right_disconnects = true
	_graph.connection_request.connect(_on_connection_request)
	_graph.disconnection_request.connect(_on_disconnection_request)
	_graph.delete_nodes_request.connect(_on_delete_nodes_request)
	_graph.popup_request.connect(_on_popup_request)
	rbox.add_child(_graph)

	_ctx_menu = PopupMenu.new()
	_ctx_menu.name = "CtxMenu"
	_ctx_menu.id_pressed.connect(_on_ctx_item_pressed)
	add_child(_ctx_menu)

	return root


func _build_rule_toolbar() -> Control:
	var vbox = VBoxContainer.new()

	# Row 1: LHS, weight, zone
	var r1 = HBoxContainer.new(); vbox.add_child(r1)
	r1.add_child(_mk("Replaces:"))
	_lhs_opt = OptionButton.new()
	_lhs_opt.custom_minimum_size = Vector2(120, 0)
	r1.add_child(_lhs_opt)
	r1.add_child(_mk("  Weight:"))
	_prob_spin = SpinBox.new()
	_prob_spin.step = 0.1; _prob_spin.min_value = 0.0; _prob_spin.max_value = 100.0
	_prob_spin.custom_minimum_size = Vector2(80, 0)
	r1.add_child(_prob_spin)
	r1.add_child(_mk("  New edge zone:"))
	_zone_opt = OptionButton.new()
	for z in ZONE_LABELS: _zone_opt.add_item(z)
	r1.add_child(_zone_opt)

	# Row 2: Condition
	var r2 = HBoxContainer.new(); vbox.add_child(r2)
	_cond_chk = CheckBox.new(); _cond_chk.text = "Condition:"
	_cond_chk.tooltip_text = "Rule only fires when this state-variable condition is true."
	_cond_chk.toggled.connect(func(on): _cond_row.visible = on)
	r2.add_child(_cond_chk)
	_cond_row = HBoxContainer.new(); _cond_row.hide(); r2.add_child(_cond_row)
	_cond_var_opt = OptionButton.new(); _cond_var_opt.custom_minimum_size = Vector2(100, 0)
	_cond_row.add_child(_cond_var_opt)
	_cond_op_opt = OptionButton.new()
	for op in OP_LABELS: _cond_op_opt.add_item(op)
	_cond_row.add_child(_cond_op_opt)
	_cond_val_sp = SpinBox.new()
	_cond_val_sp.step = 0.1; _cond_val_sp.min_value = -9999; _cond_val_sp.max_value = 9999
	_cond_val_sp.custom_minimum_size = Vector2(70, 0)
	_cond_row.add_child(_cond_val_sp)

	# Row 3: Actions
	var r3 = HBoxContainer.new(); vbox.add_child(r3)
	r3.add_child(_mk("On fire:"))
	_actions_vbox = VBoxContainer.new()
	_actions_vbox.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	r3.add_child(_actions_vbox)
	var act_add = Button.new(); act_add.text = "+ Action"
	act_add.pressed.connect(_on_add_action); r3.add_child(act_add)

	vbox.add_child(HSeparator.new())
	return vbox

# ─────────────────────────────────────────────────────────────────────────────
#  RULES — list management
# ─────────────────────────────────────────────────────────────────────────────

func _refresh_rules_list():
	_rules_list.clear()
	if not current_resource: return
	for r: GraphRule in current_resource.rules:
		_rules_list.add_item("%s → %s" % [r.lhs_symbol, r.rule_name])


func _on_add_rule():
	if not current_resource: return
	_save_current_rule()
	var rule = GraphRule.new()
	rule.rule_name = "Rule%d" % current_resource.rules.size()
	if current_resource.room_types.size() > 0:
		rule.lhs_symbol = current_resource.room_types[0].symbol
	current_resource.rules.append(rule)
	_refresh_rules_list()
	var new_idx = current_resource.rules.size() - 1
	_rules_list.select(new_idx)
	_on_rule_selected(new_idx)


func _on_delete_rule():
	if not current_resource or _rules_sel < 0: return
	current_resource.rules.remove_at(_rules_sel)
	_rules_sel = -1
	_refresh_rules_list()
	_clear_graph()


func _on_rule_selected(idx: int):
	_save_current_rule()
	_rules_sel = idx
	if not current_resource or idx < 0 or idx >= current_resource.rules.size(): return
	_load_rule_into_graph(current_resource.rules[idx])

# ─────────────────────────────────────────────────────────────────────────────
#  RULES — toolbar events
# ─────────────────────────────────────────────────────────────────────────────

func _on_add_action():
	if not current_resource or _rules_sel < 0: return
	_add_action_row("", 1.0)


func _add_action_row(var_name: String, delta: float):
	var row = HBoxContainer.new(); row.name = "ActionRow"
	var v_opt = OptionButton.new(); v_opt.name = "VarOpt"
	v_opt.custom_minimum_size = Vector2(100, 0)
	_populate_state_var_option(v_opt)
	for i in v_opt.item_count:
		if v_opt.get_item_text(i) == var_name:
			v_opt.select(i)
			break
	row.add_child(v_opt)
	var pm = Label.new(); pm.text = " ± "; row.add_child(pm)
	var d_sp = SpinBox.new(); d_sp.name = "DeltaSpin"
	d_sp.step = 0.1; d_sp.min_value = -9999; d_sp.max_value = 9999; d_sp.value = delta
	d_sp.custom_minimum_size = Vector2(70, 0); row.add_child(d_sp)
	var rm = Button.new(); rm.text = "✕"
	rm.pressed.connect(func(): row.queue_free()); row.add_child(rm)
	_actions_vbox.add_child(row)

# ─────────────────────────────────────────────────────────────────────────────
#  RULES — load / save graph ↔ rule
# ─────────────────────────────────────────────────────────────────────────────

func _load_rule_into_graph(rule: GraphRule):
	_clear_graph(); _edge_meta.clear()
	_refresh_lhs_options()
	if _lhs_opt:
		for i in _lhs_opt.item_count:
			if _lhs_opt.get_item_text(i) == rule.lhs_symbol:
				_lhs_opt.select(i)
				break
	if _prob_spin:
		_prob_spin.value = rule.probability

	var has_cond = rule.condition_var != ""
	if is_instance_valid(_cond_chk):
		_cond_chk.set_pressed_no_signal(has_cond)
	if _cond_row: _cond_row.visible = has_cond
	if is_instance_valid(_cond_var_opt):
		_populate_state_var_option(_cond_var_opt)
		for i in _cond_var_opt.item_count:
			if _cond_var_opt.get_item_text(i) == rule.condition_var:
				_cond_var_opt.select(i)
				break
	if is_instance_valid(_cond_op_opt):
		for i in _cond_op_opt.item_count:
			if _cond_op_opt.get_item_text(i) == rule.condition_op:
				_cond_op_opt.select(i)
				break
	if is_instance_valid(_cond_val_sp):
		_cond_val_sp.value = rule.condition_val

	if _actions_vbox:
		for c in _actions_vbox.get_children(): c.queue_free()
	for act in rule.actions:
		_add_action_row(act.get("var",""), act.get("delta", 1.0))

	_rebuild_ctx_menu()

	for n_data in rule.rhs_nodes:
		var gn = _make_graph_node(n_data.get("id",""), n_data.get("symbol","Symbol"))
		gn.position_offset = n_data.get("editor_pos", Vector2(100,100))
		var box = gn.get_node_or_null("Box")
		var role_row = box.get_node_or_null("RoleRow") if box else null
		var entry_chk = role_row.get_node_or_null("EntryChk") if role_row else null
		var exit_chk = role_row.get_node_or_null("ExitChk") if role_row else null
		if entry_chk: entry_chk.button_pressed = (rule.entry_node_id == n_data.get("id",""))
		if exit_chk: exit_chk.button_pressed  = (rule.exit_node_id  == n_data.get("id",""))
		_graph.add_child(gn)

	for e_data in rule.rhs_edges:
		_graph.connect_node(e_data.get("from",""), 0, e_data.get("to",""), 0)
		_edge_meta[e_data.get("from","") + "_" + e_data.get("to","")] = e_data.get("type","corridor")


func _save_current_rule():
	if not current_resource or _rules_sel < 0 or _rules_sel >= current_resource.rules.size(): return
	var rule: GraphRule = current_resource.rules[_rules_sel]

	if _lhs_opt and _lhs_opt.selected >= 0:
		rule.lhs_symbol = _lhs_opt.get_item_text(_lhs_opt.selected)
	rule.rule_name  = "Rule%d" % _rules_sel
	rule.probability = _prob_spin.value if _prob_spin else 1.0

	if is_instance_valid(_cond_chk) and _cond_chk.button_pressed and is_instance_valid(_cond_var_opt) and _cond_var_opt.selected >= 0:
		rule.condition_var = _cond_var_opt.get_item_text(_cond_var_opt.selected)
		rule.condition_op  = _cond_op_opt.get_item_text(_cond_op_opt.selected)
		rule.condition_val = _cond_val_sp.value
	else:
		rule.condition_var = ""; rule.condition_op = "<"; rule.condition_val = 0.0

	var new_actions: Array[Dictionary] = []
	if _actions_vbox:
		for row in _actions_vbox.get_children():
			var v_opt = row.get_node_or_null("VarOpt")
			var d_sp  = row.get_node_or_null("DeltaSpin")
			if v_opt and d_sp and v_opt.selected >= 0:
				new_actions.append({"var": v_opt.get_item_text(v_opt.selected), "delta": d_sp.value})
	rule.actions.clear()
	rule.actions.append_array(new_actions)

	var new_nodes: Array[Dictionary] = []
	var new_edges: Array[Dictionary] = []
	var entry_node_id = ""
	var exit_node_id = ""

	if _graph:
		for child in _graph.get_children():
			if not child is GraphNode: continue
			var sym_opt = child.get_node_or_null("Box/SymOpt")
			var symbol  = child.get_meta("symbol","Room") if not sym_opt else \
				sym_opt.get_item_text(sym_opt.selected)
			new_nodes.append({"id": child.name, "symbol": symbol, "editor_pos": child.position_offset})
			var box = child.get_node_or_null("Box")
			var role_row = box.get_node_or_null("RoleRow") if box else null
			var entry_chk = role_row.get_node_or_null("EntryChk") if role_row else null
			var exit_chk = role_row.get_node_or_null("ExitChk") if role_row else null
			if entry_chk and entry_chk.button_pressed: entry_node_id = child.name
			if exit_chk and exit_chk.button_pressed:  exit_node_id  = child.name

		for conn in _graph.get_connection_list():
			var key = str(conn.from_node) + "_" + str(conn.to_node)
			new_edges.append({"from": str(conn.from_node), "to": str(conn.to_node),
				"type": _edge_meta.get(key, "corridor")})

	rule.rhs_nodes.clear()
	rule.rhs_nodes.append_array(new_nodes)
	rule.rhs_edges.clear()
	rule.rhs_edges.append_array(new_edges)
	rule.entry_node_id = entry_node_id
	rule.exit_node_id = exit_node_id

	_refresh_rules_list()

# ─────────────────────────────────────────────────────────────────────────────
#  GRAPH NODE FACTORY
# ─────────────────────────────────────────────────────────────────────────────

func _make_graph_node(id: String, symbol: String) -> GraphNode:
	var gn = GraphNode.new()
	gn.title = symbol
	if id != "": gn.name = id
	else: _node_counter += 1; gn.name = "n%d" % _node_counter
	gn.set_slot(0, true, 0, Color.WHITE, true, 0, Color.WHITE)

	var box = VBoxContainer.new(); box.name = "Box"; gn.add_child(box)

	var sym_opt = OptionButton.new(); sym_opt.name = "SymOpt"
	_populate_palette_option(sym_opt)
	for i in sym_opt.item_count:
		if sym_opt.get_item_text(i) == symbol:
			sym_opt.select(i)
			break
	sym_opt.item_selected.connect(func(_i): _update_gn_title(gn, sym_opt))
	box.add_child(sym_opt)

	var role_row = HBoxContainer.new(); role_row.name = "RoleRow"
	var entry_chk = CheckBox.new(); entry_chk.name = "EntryChk"; entry_chk.text = "Entry"
	entry_chk.tooltip_text = "Incoming parent-graph edges connect here"
	var exit_chk = CheckBox.new(); exit_chk.name = "ExitChk"; exit_chk.text = "Exit"
	exit_chk.tooltip_text = "Outgoing parent-graph edges leave here"
	role_row.add_child(entry_chk); role_row.add_child(exit_chk)
	box.add_child(role_row)

	_tint_gn(gn, symbol)
	return gn


func _update_gn_title(gn: GraphNode, sym_opt: OptionButton):
	if sym_opt.selected < 0: return
	var sym = sym_opt.get_item_text(sym_opt.selected)
	gn.title = sym; gn.set_meta("symbol", sym); _tint_gn(gn, sym)


func _tint_gn(gn: GraphNode, symbol: String):
	if current_resource and is_instance_valid(current_resource) and current_resource.has_method("get_room_type"):
		var rt = current_resource.get_room_type(symbol)
		if rt: gn.self_modulate = rt.color.lerp(Color.WHITE, 0.5); return
	gn.self_modulate = Color.WHITE


func _refresh_graph_node_colors():
	if not _graph: return
	for child in _graph.get_children():
		if child is GraphNode:
			var sym_opt = child.get_node_or_null("Box/SymOpt")
			if sym_opt and sym_opt.selected >= 0:
				_tint_gn(child, sym_opt.get_item_text(sym_opt.selected))


func _clear_graph():
	if not _graph: return
	_graph.clear_connections()
	for c in _graph.get_children():
		if c is GraphNode: c.queue_free()

# ─────────────────────────────────────────────────────────────────────────────
#  GRAPH CALLBACKS
# ─────────────────────────────────────────────────────────────────────────────

func _on_connection_request(from_node, from_port, to_node, to_port):
	_graph.connect_node(from_node, from_port, to_node, to_port)
	_edge_meta[str(from_node)+"_"+str(to_node)] = _zone_opt.get_item_text(_zone_opt.selected)


func _on_disconnection_request(from_node, from_port, to_node, to_port):
	_graph.disconnect_node(from_node, from_port, to_node, to_port)
	_edge_meta.erase(str(from_node)+"_"+str(to_node))


func _on_delete_nodes_request(nodes):
	for nname in nodes:
		var gn = _graph.get_node_or_null(NodePath(str(nname)))
		if not gn: continue
		for conn in _graph.get_connection_list():
			if conn.from_node == nname or conn.to_node == nname:
				_graph.disconnect_node(conn.from_node, conn.from_port, conn.to_node, conn.to_port)
				_edge_meta.erase(str(conn.from_node)+"_"+str(conn.to_node))
		gn.queue_free()


func _on_popup_request(at_pos: Vector2):
	_ctx_pos = at_pos
	_ctx_menu.position = Vector2i(get_screen_position() + at_pos)
	_ctx_menu.popup()


func _on_ctx_item_pressed(id: int):
	if not current_resource: return
	var symbol = _ctx_menu.get_item_text(id)
	var gn = _make_graph_node("", symbol)
	_graph.add_child(gn)
	gn.position_offset = _graph.scroll_offset + _ctx_pos


func _rebuild_ctx_menu():
	_ctx_menu.clear()
	if not current_resource: return
	for i in current_resource.room_types.size():
		_ctx_menu.add_item(current_resource.room_types[i].symbol, i)

# ─────────────────────────────────────────────────────────────────────────────
#  STATE VARIABLES
# ─────────────────────────────────────────────────────────────────────────────

func _refresh_state_vars():
	if not _state_flow: return
	for c in _state_flow.get_children(): c.queue_free()
	if not current_resource: return
	for i in current_resource.state_variables.size():
		_add_state_var_tag(current_resource.state_variables[i], i)


func _add_state_var_tag(var_name: String, idx: int):
	var tag = HBoxContainer.new()
	var lbl = Label.new(); lbl.text = var_name; tag.add_child(lbl)
	var rm = Button.new(); rm.text = "✕"
	rm.pressed.connect(func():
		current_resource.state_variables.remove_at(idx)
		_refresh_state_vars(); _refresh_cond_and_action_options())
	tag.add_child(rm); _state_flow.add_child(tag)


func _on_add_state_var():
	if not current_resource: return
	var dlg = AcceptDialog.new(); dlg.title = "Add State Variable"
	var edit = LineEdit.new(); edit.placeholder_text = "e.g. keys"
	dlg.add_child(edit)
	dlg.confirmed.connect(func():
		var v = edit.text.strip_edges()
		if v != "" and v not in current_resource.state_variables:
			current_resource.state_variables.append(v)
			_refresh_state_vars(); _refresh_cond_and_action_options()
		dlg.queue_free())
	dlg.canceled.connect(func(): dlg.queue_free())
	add_child(dlg); dlg.popup_centered()


func _refresh_cond_and_action_options():
	if _cond_var_opt: _populate_state_var_option(_cond_var_opt)
	if _actions_vbox:
		for row in _actions_vbox.get_children():
			var v = row.get_node_or_null("VarOpt")
			if v: _populate_state_var_option(v)

# ─────────────────────────────────────────────────────────────────────────────
#  OPTION HELPERS
# ─────────────────────────────────────────────────────────────────────────────

func _refresh_lhs_options():
	if not _lhs_opt: return
	var cur = _lhs_opt.get_item_text(_lhs_opt.selected) if _lhs_opt.selected >= 0 else ""
	_populate_palette_option(_lhs_opt)
	for i in _lhs_opt.item_count:
		if _lhs_opt.get_item_text(i) == cur: _lhs_opt.select(i); return


func _populate_palette_option(opt: OptionButton):
	opt.clear()
	if not current_resource: return
	for rt: RoomType in current_resource.room_types: opt.add_item(rt.symbol)


func _populate_state_var_option(opt: OptionButton):
	var cur = opt.get_item_text(opt.selected) if opt.selected >= 0 else ""
	opt.clear()
	if not current_resource: return
	for v in current_resource.state_variables: opt.add_item(v)
	for i in opt.item_count:
		if opt.get_item_text(i) == cur: opt.select(i); return

# ─────────────────────────────────────────────────────────────────────────────
#  FULL REFRESH
# ─────────────────────────────────────────────────────────────────────────────

func _refresh_all():
	_pal_sel   = -1
	_rules_sel = -1
	_refresh_palette_list()
	_refresh_rules_list()
	_refresh_state_vars()
	_refresh_lhs_options()
	_rebuild_ctx_menu()
	_clear_graph()
	if _pal_insp_box:
		for c in _pal_insp_box.get_children(): c.queue_free()
		_pal_insp_box.add_child(_mk("Select a room type from the list to edit it."))

# ─────────────────────────────────────────────────────────────────────────────
#  UTILITY
# ─────────────────────────────────────────────────────────────────────────────

func _mk(text: String) -> Label:
	var l = Label.new(); l.text = text; return l


func _vec3i_row(val: Vector3i, setter: Callable) -> HBoxContainer:
	var row = HBoxContainer.new()
	var axes = ["x:", "y:", "z:"]; var vals = [val.x, val.y, val.z]; var spins = []
	for i in 3:
		var sp = SpinBox.new()
		sp.step = 1; sp.min_value = 1; sp.max_value = 500; sp.value = vals[i]
		sp.prefix = axes[i]; sp.custom_minimum_size = Vector2(70, 0)
		spins.append(sp); row.add_child(sp)
	for i in 3:
		spins[i].value_changed.connect(func(_v):
			setter.call(Vector3i(int(spins[0].value), int(spins[1].value), int(spins[2].value))))
	return row


func _color_icon(col: Color) -> ImageTexture:
	var img = Image.create(16, 16, false, Image.FORMAT_RGBA8)
	img.fill(col)
	return ImageTexture.create_from_image(img)
