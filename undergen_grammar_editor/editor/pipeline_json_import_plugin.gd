@tool
extends EditorImportPlugin
## Auto-import plugin for *.pipeline.json files.
## Drops a .tres resource alongside the .json so it can be assigned
## directly to UnderGenWorld3D without manual import steps.

func _get_importer_name() -> String:
	return "undergen.pipeline_json"

func _get_visible_name() -> String:
	return "UnderGen Pipeline (JSON)"

func _get_recognized_extensions() -> PackedStringArray:
	return PackedStringArray(["pipeline.json"])

func _get_save_extension() -> String:
	return "tres"

func _get_resource_type() -> String:
	return "UnderGenPipeline"

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
	# Load from JSON via C++ UnderGenPipeline
	var pipeline: Resource = UnderGenPipeline.load_from_json_file(source_file)
	if pipeline == null:
		printerr("PipelineJSONImport: Failed to parse pipeline from: ", source_file)
		return ERR_PARSE_ERROR

	# Save as .tres
	var tres_path = save_path + "." + _get_save_extension()
	var err = ResourceSaver.save(pipeline, tres_path)
	if err != OK:
		printerr("PipelineJSONImport: Failed to save .tres to ", tres_path)
		return err

	print("PipelineJSONImport: ", source_file.get_file(), " → ", tres_path.get_file())
	return OK
