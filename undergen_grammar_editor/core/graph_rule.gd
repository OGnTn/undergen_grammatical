@tool
class_name GraphRule extends Resource
## A single grammar expansion rule.
## LHS: the symbol to replace.  RHS: the mini-graph it expands into.

@export var rule_name: String = "New Rule"
@export var lhs_symbol: String = "Room"
@export var probability: float = 1.0

## Entry/exit nodes link parent-graph edges through this subgraph.
@export var entry_node_id: String = ""
@export var exit_node_id: String = ""

## Structured condition — rule only fires when true (blank var = always fire).
@export_group("Condition")
@export var condition_var: String = ""   ## State variable name (blank = no condition)
@export var condition_op: String = "<"   ## < > <= >= == !=
@export var condition_val: float = 0.0

## Actions: each entry is { "var": String, "delta": float }
@export_group("Actions")
@export var actions: Array[Dictionary] = []

## RHS graph data (serialized as plain Dictionaries).
## rhs_nodes: [{ "id", "symbol", "editor_pos" }]
## rhs_edges:  [{ "from", "to", "type" }]
@export_group("Graph Data")
@export var rhs_nodes: Array[Dictionary] = []
@export var rhs_edges: Array[Dictionary] = []
