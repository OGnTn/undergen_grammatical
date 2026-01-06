class_name LevelGrammar extends Resource

# Helper to construct the Dictionary structure for C++
func generate_graph() -> Dictionary:
	return {
		"nodes": [], 
		"edges": []
	}

# -- Helper Functions for Subclasses --

func create_node(id: String, type: String, min_size: Vector3i, max_size: Vector3i) -> Dictionary:
	return {
		"id": id,
		"type": type,
		"min_size": min_size,
		"max_size": max_size
	}

func create_edge(from_id: String, to_id: String, type: String) -> Dictionary:
	return {
		"from": from_id,
		"to": to_id,
		"type": type
	}
