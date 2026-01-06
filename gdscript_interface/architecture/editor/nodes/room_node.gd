@tool
extends GraphNode

var node_id: String

func _ready():
	title = "Room"
	# Generate unique ID
	node_id = str(randi())

func get_data() -> Dictionary:
	return {
		"id": node_id,
		"type": $TypeEdit.text, # Assuming you have this node
		"min_size": Vector3i($Dimensions/MinDimensions/WidthEdit.value, $Dimensions/MinDimensions/HeightEdit.value, $Dimensions/MinDimensions/DepthEdit.value),
		"max_size": Vector3i($Dimensions/MaxDimensions/WidthEdit.value, $Dimensions/MaxDimensions/HeightEdit.value, $Dimensions/MaxDimensions/DepthEdit.value),
		"editor_offset": position_offset
	}
