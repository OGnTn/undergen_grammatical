@tool
extends EditorPlugin

# --- Sub-editors ---
var grammar_editor_instance  # Grammar rule graph editor
var pipeline_editor_instance # PCG pipeline graph editor

# --- State ---
var _current_tab: int = 0  # 0 = Grammar, 1 = Pipeline
var _selected_world_node: Node = null

func _enter_tree():
	# Preload core scripts so global class_names resolve (avoids placeholder errors).
	preload("res://addons/undergen_grammar_editor/core/level_grammar_resource.gd")
	preload("res://addons/undergen_grammar_editor/core/room_type.gd")
	preload("res://addons/undergen_grammar_editor/core/graph_rule.gd")

	# Load Grammar editor
	grammar_editor_instance = preload(
		"res://addons/undergen_grammar_editor/editor/grammar_editor.gd"
	).new()

	# Load Pipeline editor
	pipeline_editor_instance = preload(
		"res://addons/undergen_grammar_editor/editor/pipeline_editor.gd"
	).new()

	# Wrap both in a single main container
	var main_container = _build_main_screen()
	get_editor_interface().get_editor_main_screen().add_child(main_container)
	main_container.hide()

	# Register UnderGenWorld3D as a custom node type with an icon
	add_custom_type(
		"UnderGenWorld3D",
		"Node3D",
		null,  # No GDScript override needed — it's a C++ type
		get_editor_interface().get_base_control().get_theme_icon("WorldEnvironment", "EditorIcons")
	)

func _exit_tree():
	remove_custom_type("UnderGenWorld3D")
	if grammar_editor_instance and is_instance_valid(grammar_editor_instance):
		grammar_editor_instance.queue_free()
	if pipeline_editor_instance and is_instance_valid(pipeline_editor_instance):
		pipeline_editor_instance.queue_free()

# --- Main Screen Methods ---

func _has_main_screen():
	return true

func _make_visible(visible: bool):
	var main = _get_main_container()
	if main:
		main.visible = visible

func _get_plugin_name():
	return "UnderGen"

func _get_plugin_icon():
	return get_editor_interface().get_base_control().get_theme_icon("Node3D", "EditorIcons")

# --- Object Handling (determines when the tab auto-opens) ---

func _handles(object) -> bool:
	if object is LevelGrammarResource:
		return true
	if object != null and object.get_class() == "UnderGenPipeline":
		return true
	# Auto-open when an UnderGenWorld3D node is selected in the scene
	if object != null and object.get_class() == "UnderGenWorld3D":
		return true
	return false

func _edit(object):
	if object is LevelGrammarResource:
		_current_tab = 0
		_refresh_tab()
		if grammar_editor_instance:
			grammar_editor_instance.current_resource = object

	elif object != null and object.get_class() == "UnderGenPipeline":
		_current_tab = 1
		_refresh_tab()
		if pipeline_editor_instance:
			pipeline_editor_instance.current_pipeline = object

	elif object != null and object.get_class() == "UnderGenWorld3D":
		_selected_world_node = object
		_current_tab = 1
		_refresh_tab()
		# Try to load the pipeline assigned on this node
		var pipeline = object.get("pipeline")
		if pipeline and pipeline_editor_instance:
			pipeline_editor_instance.current_pipeline = pipeline

# --- Internal Helpers ---

var _main_container: Control = null
var _tab_bar: TabBar = null
var _grammar_wrapper: Control = null
var _pipeline_wrapper: Control = null

func _build_main_screen() -> Control:
	_main_container = VBoxContainer.new()
	_main_container.name = "UnderGenMainScreen"
	_main_container.set_anchors_preset(Control.PRESET_FULL_RECT)
	_main_container.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	_main_container.size_flags_vertical = Control.SIZE_EXPAND_FILL

	# --- Tab Bar ---
	_tab_bar = TabBar.new()
	_tab_bar.add_tab("Level Grammar")
	_tab_bar.add_tab("Physical Pipeline")
	_tab_bar.tab_changed.connect(_on_tab_changed)
	_main_container.add_child(_tab_bar)

	# --- Grammar Editor Wrapper ---
	_grammar_wrapper = grammar_editor_instance
	grammar_editor_instance.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	grammar_editor_instance.size_flags_vertical = Control.SIZE_EXPAND_FILL
	_main_container.add_child(grammar_editor_instance)

	# --- Pipeline Editor Wrapper ---
	_pipeline_wrapper = pipeline_editor_instance
	pipeline_editor_instance.size_flags_horizontal = Control.SIZE_EXPAND_FILL
	pipeline_editor_instance.size_flags_vertical = Control.SIZE_EXPAND_FILL
	pipeline_editor_instance.hide()
	_main_container.add_child(pipeline_editor_instance)

	return _main_container

func _get_main_container() -> Control:
	return _main_container

func _on_tab_changed(tab_idx: int):
	_current_tab = tab_idx
	_refresh_tab()

func _refresh_tab():
	if grammar_editor_instance:
		grammar_editor_instance.visible = (_current_tab == 0)
	if pipeline_editor_instance:
		pipeline_editor_instance.visible = (_current_tab == 1)
	if _tab_bar:
		_tab_bar.current_tab = _current_tab
