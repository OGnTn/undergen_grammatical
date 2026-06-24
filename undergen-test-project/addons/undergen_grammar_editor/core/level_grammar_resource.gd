@tool
class_name LevelGrammarResource extends Resource
## Top-level grammar resource.
## Holds the room type palette, state variable names, and expansion rules.

@export var axiom: String = "Start"
@export var room_types: Array[RoomType] = []
@export var state_variables: Array[String] = []
@export var rules: Array[GraphRule] = []

# --- Palette helpers ---

func get_room_type(symbol: String) -> RoomType:
	for rt in room_types:
		if rt.symbol == symbol:
			return rt
	return null

func get_palette_symbols() -> Array[String]:
	var out: Array[String] = []
	for rt in room_types:
		out.append(rt.symbol)
	return out

func get_palette_color(symbol: String) -> Color:
	var rt = get_room_type(symbol)
	return rt.color if rt else Color(0.5, 0.5, 0.5)

# --- Rule helpers ---

func get_rules_for_symbol(symbol: String) -> Array[GraphRule]:
	var matching: Array[GraphRule] = []
	for r in rules:
		if r.lhs_symbol == symbol:
			matching.append(r)
	return matching
