@tool
extends EditorImportPlugin
## Auto-import plugin for *.grammar.json files.
## Drops a .tres resource alongside the .json so it can be used directly
## by UnderGenGrammarNode without manual import steps.

func _get_importer_name() -> String:
	return "undergen.grammar_json"

func _get_visible_name() -> String:
	return "UnderGen Grammar (JSON)"

func _get_recognized_extensions() -> PackedStringArray:
	return PackedStringArray(["grammar.json"])

func _get_save_extension() -> String:
	return "tres"

func _get_resource_type() -> String:
	return "LevelGrammarSpec"

func _get_preset_count() -> int:
	return 0

func _get_import_order() -> int:
	return 0

func _get_priority() -> float:
	return 1.0

func _get_option_visibility(path: String, option_name: StringName, options: Dictionary) -> bool:
	return true

func _import(source_file: String, save_path: String, options: Dictionary,
			 platform_variants: Array[String], gen_files: Array[String]) -> Error:
	# Read the JSON file
	var file = FileAccess.open(source_file, FileAccess.READ)
	if file == null:
		printerr("GrammarJSONImport: Cannot open source file: ", source_file)
		return ERR_FILE_CANT_OPEN
	var text = file.get_as_text()
	file.close()

	# Parse via C++ LevelGrammarSpec (handles JSON → typed properties correctly)
	var spec: RefCounted = ClassDB.instantiate("LevelGrammarSpec")
	if spec == null:
		printerr("GrammarJSONImport: LevelGrammarSpec not available — is the GDExtension loaded?")
		return ERR_UNAVAILABLE

	var err = spec.from_json(text)
	if err != OK:
		printerr("GrammarJSONImport: JSON parse error in ", source_file)
		return err

	# Save as .tres
	var tres_path = save_path + "." + _get_save_extension()
	err = ResourceSaver.save(spec, tres_path)
	if err != OK:
		printerr("GrammarJSONImport: Failed to save .tres to ", tres_path)
		return err

	print("GrammarJSONImport: ", source_file.get_file(), " → ", tres_path.get_file())
	return OK
