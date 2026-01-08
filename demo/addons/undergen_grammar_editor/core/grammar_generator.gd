class_name GrammarGenerator extends RefCounted

var grammar: LevelGrammarResource
var rng: RandomNumberGenerator

func _init(p_grammar: LevelGrammarResource, seed_val: int):
    grammar = p_grammar
    rng = RandomNumberGenerator.new()
    rng.seed = seed_val

# Generates the final node/edge lists for LevelDensityGrid
func generate(iterations: int = 4, initial_state: Dictionary = {}, max_nodes: int = 100) -> Dictionary:
    # 1. Initialize State & Axiom
    var state = initial_state.duplicate()
    var nodes = []
    var edges = []
    
    var root_node = {
        "id": "root",
        "symbol": grammar.axiom,
        "type": "generic", # Mapping symbol->type?
        "min_size": Vector3i(5, 5, 5),
        "max_size": Vector3i(10, 10, 10),
        "constraints": {"fixed_pos": Vector3(0, 0, 0)}
    }
    nodes.append(root_node)
    
    # 2. Iterative Expansion
    for i in range(iterations):
        if nodes.size() >= max_nodes:
            print("Grammar: Reached max node count (%d). Stopping expansion." % max_nodes)
            break
            
        var replacements = {} # old_id -> { "entry": new_id, "nodes": [], "edges": [] }
        var preserved_nodes = []
        var any_change = false
        
        # A. Determine replacements
        for node in nodes:
            var symbol = node["symbol"]
            var candidates = grammar.get_rules_for_symbol(symbol)
            var rule = _pick_rule(candidates, state)
            
            if rule:
                any_change = true
                
                # Apply Actions to State
                for k in rule.actions:
                    var val = rule.actions[k]
                    if val is int or val is float:
                        state[k] = state.get(k, 0) + val
                    else:
                        state[k] = val
                        
                # Instantiate Subgraph
                var sub_nodes = []
                var sub_edges = []
                var id_map = {} # local_rule_id -> global_unique_id
                
                # Create Nodes
                for rhs_n in rule.rhs_nodes:
                    var unique_id = node["id"] + "-" + rhs_n["id"] + "-" + str(i)
                    id_map[rhs_n["id"]] = unique_id
                    
                    var new_n = {
                        "id": unique_id,
                        "symbol": rhs_n["symbol"],
                        "type": rhs_n.get("symbol", "generic").to_lower(),
                        "constraints": rhs_n.get("constraints", {}).duplicate()
                    }
                    if rhs_n.has("min_size"): new_n["min_size"] = rhs_n["min_size"]
                    if rhs_n.has("max_size"): new_n["max_size"] = rhs_n["max_size"]
                    if rhs_n.has("vox_path"): new_n["vox_path"] = rhs_n["vox_path"]
                    if rhs_n.has("connect_to"): new_n["connect_to"] = rhs_n["connect_to"]
                    
                    # Resolve Relative Constraints
                    var c = new_n["constraints"]
                    if c.has("relative_to"):
                        var rel = c["relative_to"]
                        if rel == "LHS" or rel == "SELF":
                            # Relative to the node we are replacing (usually invalid as it's gone? 
                            # Or relative to its *original* position? No, grammar is topological)
                            # Actually, usually relative to "PREV" or "PARENT". 
                            # If replacing A -> B, 'B' is new. A is gone. 
                            # If we want B relative to A's *parent*, that's hard.
                            # But usually we define relativity *internal* to the rule.
                            # e.g. Rule: A -> B + C. C relative to B.
                            pass
                        elif id_map.has(rel):
                            # Internal reference (C relative to B)
                            c["relative_to"] = id_map[rel]
                            
                    sub_nodes.append(new_n)
                
                # Create Internal Edges
                for rhs_e in rule.rhs_edges:
                    sub_edges.append({
                        "from": id_map[rhs_e["from"]],
                        "to": id_map[rhs_e["to"]],
                        "type": rhs_e.get("type", "corridor") # Propagate Type
                    })
                
                # Determine Entry/Exit Points
                var entry_id = sub_nodes[0]["id"] if sub_nodes.size() > 0 else node["id"]
                var exit_id = entry_id # Default to same node if not specified
                
                if rule.entry_node_id != "" and id_map.has(rule.entry_node_id):
                    entry_id = id_map[rule.entry_node_id]
                    
                if rule.exit_node_id != "" and id_map.has(rule.exit_node_id):
                    exit_id = id_map[rule.exit_node_id]

                
                replacements[node["id"]] = {
                    "entry": entry_id,
                    "exit": exit_id,
                    "nodes": sub_nodes,
                    "edges": sub_edges
                }
            else:
                preserved_nodes.append(node)
                
        if not any_change:
            break
            
        # B. Rebuild Graph
        var next_nodes = []
        var next_edges = []
        
        # Add preserved nodes
        next_nodes.append_array(preserved_nodes)
        
        # Add replacement nodes & internal edges
        for key in replacements:
            var data = replacements[key]
            next_nodes.append_array(data["nodes"])
            next_edges.append_array(data["edges"])
            
        # Rewire old edges
        # If an edge connected NodeA -> NodeB:
        # Case 1: Both preserved. Keep edge.
        # Case 2: NodeA replaced. Edge comes from NodeA_Entry.
        # Case 3: NodeB replaced. Edge goes to NodeB_Entry.
        
        for old_edge in edges:
            var from = old_edge["from"]
            var to = old_edge["to"]
            
            var new_from = from
            if replacements.has(from):
                new_from = replacements[from]["exit"]
            
            var new_to = to
            if replacements.has(to):
                new_to = replacements[to]["entry"]
                
            if new_from != new_to: # Prevent self-loops if degenerate
                next_edges.append({
                    "from": new_from,
                    "to": new_to,
                    "type": old_edge["type"]
                })
                
        nodes = next_nodes
        edges = next_edges
        
    # 3. Post-Process: Resolve Connections (LOOPS)
    print("Grammar: Resolving Connections...")
    var connection_edges = []
    
    for i in range(nodes.size()):
        var node = nodes[i]
        if node.has("connect_to"):
            var ct = node["connect_to"]
            var target_sym = ct.get("target", "")
            
            var best_idx = -1
            var candidates = []
            
            # Find Candidates
            for j in range(nodes.size()):
                if i == j: continue
                var other = nodes[j]
                
                # Check symbol match (case insensitive?)
                if other.get("symbol", "") == target_sym or target_sym == "":
                    # Check if edge already exists? (Expensive O(E))
                    # For now assume multigraph is OK or handled by logic
                    # Calculate Logic Distance (Topology Distance? No, physics distance is unknown yet)
                    # We can only reliably use "Random" or "First" at this stage really..
                    # UNLESS user means "Closest in Graph Topology" (BFS)
                    # OR we wait until C++ physics step to resolve this?
                    # BUT we need to output EDGES for the solver.
                    # Wait. If we are in Grammar Generator, we have NO spatial information.
                    # We only create logical edges.
                    # "Closest" implies spatial.
                    # This logic FLIPS functionality.
                    # If the user wants a loop, they want a LOGICAL connection.
                    # So "Closest" is ambiguous.
                    # Alternative: "Random" valid target.
                    candidates.append(j)
            
            if candidates.size() > 0:
                best_idx = candidates[rng.randi() % candidates.size()]
                
                if best_idx != -1:
                    print("Grammar: Connected %s -> %s (Loop)" % [node["id"], nodes[best_idx]["id"]])
                    connection_edges.append({
                        "from": node["id"],
                        "to": nodes[best_idx]["id"],
                        "type": "corridor" # Default type
                    })

    edges.append_array(connection_edges)

    return {"nodes": nodes, "edges": edges}

func _pick_rule(candidates: Array, state: Dictionary) -> GraphRule:
    var valid_rules = []
    
    # Filter by Condition
    for r in candidates:
        if r.condition == "":
            valid_rules.append(r)
        else:
            if _evaluate_condition(r.condition, state):
                valid_rules.append(r)
    
    if valid_rules.is_empty():
        return null
        
    var roll = rng.randf()
    var accum = 0.0
    # Normalize probabilities? Or just use raw check.
    # Let's sum valid probs first to normalize
    var total_prob = 0.0
    for r in valid_rules: total_prob += r.probability
    
    if total_prob == 0: return valid_rules[0]
    
    var pick = rng.randf() * total_prob
    for r in valid_rules:
        accum += r.probability
        if pick <= accum:
            return r
    return valid_rules[0]

func _evaluate_condition(cond: String, state: Dictionary) -> bool:
    # Simple expression parser: "key < 3"
    # Supported ops: <, >, <=, >=, ==, !=
    # Format: "VAR OP VALUE"
    var parts = cond.split(" ")
    if parts.size() != 3:
        print("Grammar: Invalid condition format '%s'" % cond)
        return false
        
    var key = parts[0]
    var op = parts[1]
    var val_str = parts[2]
    var val = val_str.to_float()
    
    var current_val = state.get(key, 0)
    
    match op:
        "<": return current_val < val
        ">": return current_val > val
        "<=": return current_val <= val
        ">=": return current_val >= val
        "==": return current_val == val
        "!=": return current_val != val
        
    return false
