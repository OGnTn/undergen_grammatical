@tool
extends EditorPlugin

var editor_instance

func _enter_tree():
    # Load the editor scene (or script class)
    editor_instance = preload("res://addons/undergen_grammar_editor/editor/grammar_editor.gd").new()
    # Add to the editor viewport (Dock or Main Screen)
    # Main Screen is better for GraphEdit
    get_editor_interface().get_editor_main_screen().add_child(editor_instance)
    editor_instance.hide()

func _exit_tree():
    if editor_instance:
        editor_instance.queue_free()

func _has_main_screen():
    return true

func _make_visible(visible):
    if editor_instance:
        editor_instance.visible = visible

func _get_plugin_name():
    return "Level Grammar"

func _get_plugin_icon():
    return get_editor_interface().get_base_control().get_theme_icon("GraphEdit", "EditorIcons")

func _handles(object):
    # Robust check in case class_name cache is stale
    if object is LevelGrammarResource:
        return true
    # Fallback check
    if object.get_script() and object.get_script().resource_path.ends_with("level_grammar_resource.gd"):
        return true
    return false

func _edit(object):
    if editor_instance:
        editor_instance.current_resource = object
