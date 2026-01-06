@tool
class_name LevelGraphResource extends Resource

# We save the raw dictionary data (Nodes + Edges) here
# so the C++ generator can read it directly.
@export var graph_data: Dictionary = {}

# We also save the editor-only data (node positions, offsets)
# so you can re-open the graph and see it nicely layouted.
@export var editor_metadata: Dictionary = {}
