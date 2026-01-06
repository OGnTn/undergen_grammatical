class_name LevelGrammarResource extends Resource

@export var axiom: String = "Start"
@export var rules: Array[GraphRule] = []

func get_rules_for_symbol(symbol: String) -> Array[GraphRule]:
    var matching: Array[GraphRule] = []
    for r in rules:
        if r.lhs_symbol == symbol:
            matching.append(r)
    return matching
