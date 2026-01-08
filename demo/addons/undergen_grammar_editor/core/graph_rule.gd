class_name GraphRule extends Resource

@export var rule_name: String = "New Rule"
@export var lhs_symbol: String = "Room" # The symbol to replace
@export var probability: float = 1.0
@export var entry_node_id: String = ""
@export var exit_node_id: String = ""

# State Logic
@export var condition: String = "" # e.g. "keys < 3"
@export var actions: Dictionary = {} # e.g. {"keys": 1} (Adds 1)

# The RHS (Right Hand Side) is a mini-graph
# Instead of storing raw dictionaries, we can store a list of Node definitions and Connections.
# We will use simple Dictionaries for storage to keep it serializable easily, or custom internal classes.
@export var rhs_nodes: Array[Dictionary] = [] # { "id": "local_1", "symbol": "Type", "constraints": {} }
@export var rhs_edges: Array[Dictionary] = [] # { "from": "local_1", "to": "local_2" }

# Helper to add a node to RHS
func add_node(local_id: String, symbol: String, constraints: Dictionary = {}) -> void:
    rhs_nodes.append({
        "id": local_id,
        "symbol": symbol,
        "constraints": constraints
    })

func add_edge(from_id: String, to_id: String) -> void:
    rhs_edges.append({
        "from": from_id,
        "to": to_id
    })
