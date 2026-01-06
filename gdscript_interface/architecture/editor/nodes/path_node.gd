@tool
extends GraphNode

func get_data() -> Dictionary:
	return {
		"type": $TypeEdit.text
	}
