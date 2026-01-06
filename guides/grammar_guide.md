# Level Grammar System: Conceptual Guide

## 1. Nodes = Volumes (Rooms)
In this system, every **Node** in the graph represents a **Volume of Space** in the final level. 
*   We call them **"Room Nodes"** in the editor because the layout solver eventually turns them into physical rooms (or caves, or junctions).
*   **Abstract Symbols**: Initially, a node might be abstract (e.g., `BossArena`).
*   **Concrete Rooms**: If no more rules apply to a `BossArena` node (it is a **Leaf** or **Terminal**), the generator creates a physical Room of that type.

## 2. Edges = Connections (Implicit Paths)
The **Edges** (white lines) represent topological connections. The solver ensures these nodes are reachable from each other.
*   **Doorways**: If nodes are adjacent, the edge becomes a door.
*   **Carved Paths**: If nodes are far apart, the edge becomes a simple carved path.

### Note on "Corridor Nodes"
Sometimes you want a **specific room** that acts as a hallway (e.g., a "Grand Hall" or "Trapped Corridor"). In that case, you add a **Node** for it.
*   `[Room A] -- [Room B]` = A path directly connects them.
*   `[Room A] -- [CorridorNode] -- [Room B]` = Room A connects to a Corridor Room, which connects to Room B.

### Path Zones
Even implicit paths are smart! They are strictly marked with the Zone ID **"corridor"** (or whatever type you set on the edge).
*   You can set your Roaming Enemy `SpawnableObject` to require `"corridor"`.
*   The `Populator` will automatically restict them to hallways.

## 3. Implicit "Leaf" Nodes
There is no explicit "Leaf Node" button. It is determined by your Rules.
*   **Non-Terminal**: A symbol is "Non-Terminal" if you have a Rule for it (e.g., `LHS: Dungeon`). The system will replace it.
*   **Terminal (Leaf)**: A symbol is "Terminal" if **no rules match it** (e.g., `LHS: TreasureRoom` has no rule). The system stops there and builds it.

## 4. Example Workflow
1.  **Rule 1**: `Start` -> `[EntryRoom] -- [Hallway] -- [BossAntechamber]`
    *   `Start` is replaced.
    *   `EntryRoom`, `Hallway`, `BossAntechamber` are created.
2.  **Rule 2**: `Hallway` -> `[CorridorSegment] -- [EnemyAmbush] -- [CorridorSegment]`
    *   The `Hallway` is replaced by a longer chain.
3.  **Result**: `EntryRoom`, `CorridorSegment`, `EnemyAmbush`, `BossAntechamber` are likely **Leaves** (no more rules). They become the final level geometry.

## 5. State & Logic
The grammar is "Stateful", meaning it can remember variables (counters, flags) and use them to decide which rules to run.

### Probability
*   **What it is**: The chance this rule is picked if multiple rules match the same symbol.
*   **Usage**: If you have two rules for `Room`, one with `Prob: 1.0` and one with `Prob: 0.5`, the standard room is twice as likely.

### Condition
*   **What it is**: A requirement for the rule to run.
*   **Format**: `VARIABLE OPERATOR VALUE`. (e.g., `keys >= 3`, `boss_dead == 1`).
*   **Usage**: Create a "BossRoom" rule that only runs if `keys >= 3`.

### Action Effects
*   **What it is**: Changes to state variables when the rule is applied.
*   **Format**: JSON Dictionary. `{"key_name": 1}` adds 1 to a number. `{"flag": 1}` sets a value.
*   **Usage**: A "KeyRoom" rule might have Action `{"keys": 1}`. When the generator places this room, it increments the global key counter, potentially unlocking the Boss Room rule later in generation.

## 6. Common Recipes

### The "Key & Chest" Pattern
You want the generator to place a Key Room *before* it tries to place a Chest Room (which requires a key).

**1. The Key Room (Giver)**
*   **LHS**: `KeyLocation`
*   **RHS**: `[KeyRoomNode]` (Terminal)
*   **Action**: `{"keys": 1}`
*   *Effect*: When this rule runs, it places a room and increments the `keys` counter.

**2. The Chest Room (Taker)**
*   **LHS**: `TreasureLocation`
*   **RHS**: `[ChestRoomNode]` (Terminal)
*   **Condition**: `keys >= 1`
*   **Action**: `{"keys": -1}` (Optional, if the key is consumed)
*   *Effect*: This rule is **only valid** if a Key Room has already been generated.

**3. The Workflow**
*   Your main graph expands. It places some `KeyLocation` symbols.
*   The generator picks the `KeyLocation` rule -> `keys` becomes 1.
*   Later, it encounters a `TreasureLocation`.
*   It checks rules. The `ChestRoom` rule sees `keys >= 1` (True!) and runs.
*   *Result*: A guaranteed valid key-chest pair.

## 7. Troubleshooting: Infinite Recursion
If your level has way too many rooms, check for self-referencing rules!

**The Recursive Mistake:**
*   Rule: `Start` -> `[Start] -- [End]`
*   *Iteration 1*: `Start` + `End`
*   *Iteration 2*: (`Start` + `End`) + `End`
*   *Iteration 3*: ((`Start` + `End`) + `End`) + `End`
*   *Result*: The `Start` node is never "finished" because it keeps replacing itself with another copy of `Start`.

**The Solution: Terminals**
Use a different symbol for the result!
*   Rule: `Start` -> `[EntryRoom] -- [EndRoom]`
*   Ensure `EntryRoom` has **no rules** (Terminal).
*   Ensure `EndRoom` has **no rules** (Terminal).
*   *Iteration 1*: `Entry` + `End`
*   *Iteration 2*: Nothing happens. `Entry` matches no rules. `End` matches no rules. Generation stops perfectly.
