treat the level as a **semantic composition of spatial forms**, not as a simulation of what the place is used for.

The relevant question is not “where should the extraction room go?” but:

> What spatial role does this volume play in the player’s experience and in the composition of the surrounding geometry?

## A vocabulary of shape semantics

Instead of generic “rooms,” define volumes by geometric function.

### Void forms

* **Primary void** — the dominant open volume organizing the level
* **Pocket** — a small enclosed volume attached to another space
* **Bay** — a lateral expansion that remains visually connected to the main space
* **Throat** — a compressed connector between larger spaces
* **Shaft** — strongly vertical void
* **Slot** — narrow void with tall opposing walls
* **Bowl** — space widening upward from a lower center
* **Dome** — enclosed volume with a prominent ceiling
* **Chasm** — negative space primarily experienced from its boundaries
* **Undercroft** — navigable or visible space beneath another route
* **Gallery** — elongated space with a consistent lateral relationship
* **Reveal chamber** — space whose purpose is to expose a larger composition

### Positive-mass forms

The image is defined just as much by what was **not carved**.

* **Anchor mass** — dominant preserved rock body
* **Occluder** — mass blocking a complete view and producing partial reveals
* **Divider** — separates nearby spaces without fully disconnecting them
* **Buttress** — wall mass projecting into a void
* **Spine** — elongated mass separating parallel paths
* **Island** — isolated mass inside a larger void
* **Canopy** — overhead mass compressing part of an otherwise large chamber
* **Screen** — thin or perforated mass creating layered silhouettes

The large central pillar in the reference is simultaneously an anchor, occluder, divider, and navigation landmark.

That is the kind of multi-role geometry worth generating.

---

## Describe a level as a spatial sentence

The reference could be abstracted as:

> A compressed entrance reveals an asymmetric primary void wrapped around a large occluding anchor mass. Traversal occupies several elevation bands along the void’s boundary. Secondary bays penetrate the surrounding shell. Cross-void connections are rare, while vertical and diagonal sightlines are common.

That description contains no mine semantics, but it captures much of the scene.

A generator could represent this as:

```text
PrimaryVoid:
    shape: elongated_asymmetric
    verticality: high
    openness: 0.8

AnchorMass:
    relation: inside PrimaryVoid
    eccentricity: 0.25
    occlusion_target: 0.45
    route_wrap: 180–270 degrees

ElevationBands:
    count: 4
    spacing: irregular
    overlap_in_view: high

SecondaryVoids:
    arrangement: clustered
    attachment: boundary
    depth_variation: high

Crossings:
    frequency: low
    prominence: high
```

This becomes a specification for the density field rather than a theme specification.

---

# Separate topology, embedding, and geometry

These are three different generation problems and should remain separate.

## 1. Topology

Determine abstract connectivity:

* Which spaces connect?
* Which paths form loops?
* Which routes rejoin?
* Which areas are optional?
* Which spaces overlook each other?
* Which routes pass above or below others?
* Which locations are visible but not directly reachable?

At this stage, nodes do not need precise coordinates or dimensions.

You might use nodes such as:

```text
Reveal
Hub
Pocket
Overlook
Crossing
VerticalConnector
DeadEnd
LoopReturn
VistaTarget
```

These are experiential and geometric roles, not thematic rooms.

Edges can carry relational properties:

```text
Edge:
    traversal_type: boundary | interior | crossing | vertical
    enclosure: 0..1
    width_profile: narrow → wide
    curvature: low | medium | high
    exposure: 0..1
    destination_visibility: hidden | partial | continuous
```

## 2. Spatial embedding

Place the graph in 3D while satisfying relationships:

* Overlooks must face spaces below.
* Parallel routes should occasionally become mutually visible.
* A reveal node should sit behind an occluder or throat.
* A crossing should span meaningful negative space.
* A loop return should reconnect from a noticeably different elevation or direction.
* Important spaces should not all occupy the same plane.

## 3. Geometric realization

Only after the embedding is stable do you turn nodes and edges into SDF or density operations:

* Void primitives
* Preserved masses
* Path sweeps
* Ledge cuts
* Ceiling expansions
* Undercuts
* Noise
* Warping

This prevents the density field from being responsible for deciding the level’s composition.

---

# Generate relationships, not isolated rooms

A chamber is rarely interesting by itself. Its relationship to other spaces is what makes it useful.

Useful spatial relationships include:

### Containment

A small chamber sits within the perceived envelope of a much larger cavern.

### Adjacency

Two spaces share a boundary but differ in floor height, width, or enclosure.

### Overlap

Spaces occupy overlapping horizontal coordinates at different elevations.

### Intervisibility

Two spaces can see each other but have no immediate connection.

### Occluded adjacency

Spaces are physically close but separated visually by a mass.

### Suspension

A route or space sits inside a void without touching its floor.

### Wrap

A path follows part of the boundary of a mass or void.

### Penetration

A path passes through a mass that otherwise acts as a divider.

### Bridging

A path crosses open negative space rather than following its boundary.

### Convergence

Several routes become visible together before physically meeting.

The image uses wrap, overlap, suspension, intervisibility, and occlusion heavily.

You can represent these as explicit constraints in your spatial graph:

```text
Route A wraps AnchorMass
Route B overlooks PrimaryVoid
Route C is visible from A but unreachable from A
Pocket D is occluded until the final 8 metres of its approach
Bridge E crosses PrimaryVoid at its narrowest upper section
```

---

# Use macro-geometry before path carving

A common procedural approach is:

1. Place rooms.
2. Connect them.
3. Add noise.

For this type of space, invert part of that process:

1. Define the dominant void.
2. Define the dominant preserved masses.
3. Establish the boundary regions available for traversal.
4. Embed paths and secondary volumes relative to those forms.
5. Locally modify the macro-geometry to support them.

The reference’s central rock pillar should not emerge accidentally from room placement. It should be an intentional **positive-space primitive**.

A useful macro construction could be:

```text
world solid
- primary cavern void
+ central anchor mass
- secondary wall bays
- vertical slot
- lower chasm
- upper ceiling pockets
- route clearances
```

The ordering matters. The route geometry should adapt to the macro-form rather than erase it.

---

# Represent spatial character as continuous fields

You can place low-resolution semantic fields over the level independently from density.

## Openness field

How much empty space should surround a point?

* High values expand caverns.
* Low values compress tunnels.
* Gradients create natural transitions.

## Verticality field

How strongly should the geometry emphasize the vertical axis?

* High verticality produces shafts, tall walls, stacked routes, and narrow horizontal extents.
* Low verticality produces broad chambers and lateral galleries.

## Exposure field

How much should a route border open negative space?

* High exposure creates cliff paths, balconies, suspended sections, and distant views.
* Low exposure creates enclosed corridors.

## Enclosure field

Related to exposure but not identical. A route can be open on one side while strongly enclosed above and behind.

## Occlusion field

Where should preserved mass interrupt long sightlines?

This can encourage:

* Central pillars
* Wall projections
* Ceiling drops
* Route bends
* Partial reveals

## Prominence field

Where should geometry produce memorable silhouettes?

Use it to reduce clutter, increase local scale contrast, or preserve characteristic masses.

## Connectivity pressure

Where is the generator encouraged to create secondary connections, loops, or openings?

## Surface suitability

Where can traversable floors or ledges plausibly exist without destroying the macro-form?

These fields can guide your carving parameters without being part of the final density function.

---

# Give paths changing spatial states

Do not treat an entire Bezier connection as one uniform tunnel. Define a sequence of spatial states along it.

For example:

```text
enclosed
→ compressed
→ partial side opening
→ exposed ledge
→ widened overlook
→ bridge
→ enclosed arrival
```

Each state controls the cross-section and its relationship with existing voids.

A path sample could contain:

```text
PathState:
    width
    height
    floor_flatness
    ceiling_offset
    left_wall_presence
    right_wall_presence
    exposure
    lateral_bias
    bank
    local_noise_scale
    connection_to_existing_void
```

Then interpolate between state keys along the spline.

### Enclosed state

Carve a mostly complete cross-section.

### Boundary state

Carve into an existing void boundary, preserving one side as open.

### Ledge state

Carve a floor and inner clearance but avoid carving the outer side.

### Balcony state

Expand locally toward an existing void and flatten the floor.

### Bridge state

Carve only endpoint clearances; let a separate surface generator create the span.

### Crevice state

Use a tall, narrow anisotropic profile.

This creates meaningful variation without requiring more random noise.

---

# Use contextual morphological operators

Think of generation as applying verbs to existing geometry.

## Inflate

Expand a void locally.

Useful for arrival spaces, overlooks, or visual breathing room.

## Pinch

Reduce width or height.

Useful before reveals or transitions.

## Split

Insert a preserved mass through a void, creating two partial channels.

## Wrap

Bias a path around a positive mass or along a void boundary.

## Terrace

Create multiple gravity-aligned ledges at related heights.

## Undercut

Remove material beneath a shelf or route.

This increases exposure and produces stronger silhouettes.

## Perforate

Cut one or more partial openings through a dividing mass.

## Occlude

Add or preserve mass between a viewpoint and destination.

## Frame

Shape nearby mass so that a distant feature occupies a controlled opening.

## Cascade

Place related spaces at descending elevations with overlapping footprints.

## Braid

Allow two routes to approach, diverge, cross in view, and reconnect.

## Stack

Place routes or spaces at several elevations around a common void.

These operators can be applied based on graph relationships and local fields rather than randomly.

---

# Generate around void boundaries

The scene’s paths largely occupy the **interface between solid rock and major empty space**.

After creating a primary void, extract a coarse representation of its boundary. From that boundary you can identify candidate regions for:

* Ledges
* Shelves
* Bays
* Overlooks
* Wall paths
* Bridges between nearby boundary sections
* Vertical stacks of paths
* Recesses behind buttresses

For each boundary point or patch, calculate:

```text
surface normal
distance to void center
height above void floor
local curvature
overhang
visible open-space solid angle
distance to opposing wall
distance to existing route
```

These measurements allow spatially meaningful decisions.

For example:

```text
high open-space angle
+ approximately vertical wall
+ moderate opposing-wall distance
= good overlook or ledge candidate

two boundary patches facing one another
+ compatible elevation
+ substantial empty space between them
= bridge candidate

concave boundary patch
+ low exposure
= pocket or sheltered route candidate

convex buttress
+ open space on both sides
= wraparound route or reveal candidate
```

---

# Explicitly design the positive-to-negative-space ratio

A cave generator often focuses almost entirely on empty space. But scenes like this depend on the arrangement of remaining solid masses.

You can measure properties such as:

* Percentage of the primary void occupied by internal masses
* Number of major isolated silhouettes
* Occlusion percentage from key viewpoints
* Average angular size of the dominant mass
* Ratio of boundary-following path to free-spanning path
* Number of visible elevation layers
* Distribution of empty space above and below traversal

For the reference, you might aim for:

```text
one dominant internal mass
two to four secondary buttresses
high vertical empty volume
moderate visual occlusion
high path overlap in screen space
low number of full-width crossings
high boundary-following route percentage
```

Without the central positive mass, the same path layout would read as a generic open cavern.

---

# Use visibility as part of shape generation

Visibility should influence the density field, not merely be evaluated afterward.

At important points, cast a coarse set of rays and estimate:

* Visible empty volume
* Percentage of destination visible
* Number of route layers visible
* Dominant silhouette complexity
* Depth range
* Degree of enclosure
* Sky or ceiling visibility
* Amount of foreground occlusion

Then adjust geometry.

Examples:

```text
Reveal point sees too much before arrival:
    thicken occluder
    narrow throat
    bend approach
    lower ceiling

Overlook feels visually flat:
    undercut foreground ledge
    expose lower routes
    move opposing wall farther away
    introduce an intermediate silhouette

Primary void lacks depth:
    add overlapping buttress
    stagger secondary pockets
    vary route distances from camera
```

This is essentially procedural composition.

---

# Generate hierarchical scale contrast

The reference combines several scales:

1. **Macro:** primary cavern and central pillar
2. **Meso:** bays, terraces, shafts, ledges
3. **Route:** paths, stairs, bridges
4. **Micro:** rock breakup and surface noise

Each scale should have its own generation pass and frequency range.

Do not expect warped noise to create the meso-scale forms. Noise can roughen a wall, but it rarely creates a convincing terrace, buttress, reveal, or layered chasm.

A useful ordering is:

```text
macro voids and masses
→ meso spatial operators
→ route carving
→ local adaptation
→ geological displacement
→ material detail
```

Keep the amplitudes separated. Micro-noise should never erase the intent of the macro or meso forms.

---

# A shape-oriented generator for this reference

One possible abstract recipe:

## Macro composition

* Create one tall asymmetric primary void.
* Place one large off-centre anchor mass inside it.
* Preserve a deep lower negative-space region.
* Add a heavier ceiling mass on one side to make the cavern asymmetric.
* Cut three clustered bays into the surrounding shell.

## Traversal composition

* Generate three main elevation bands.
* Wrap the middle band around 50–70% of the anchor mass.
* Let the upper band cross behind the anchor.
* Let the lower band appear intermittently beneath the primary path.
* Create only one major cross-void connection.
* Add several visible but disconnected ledges.

## Reveal structure

* Begin in a compressed boundary passage.
* Reveal only part of the anchor mass.
* Later reveal the lower chasm.
* Later still expose the upper route and opposite-side bays.
* Ensure no single viewpoint explains the entire level.

## Shape variation

* Widen paths where they face the primary void.
* Compress them behind the anchor.
* Increase ceiling height near crossings.
* Create undercuts beneath prominent overlooks.
* Add vertical slots where several elevation bands overlap.
* Preserve smooth navigable floors but warp surrounding walls according to a larger geological frame.

---

# A compact data model

You could encode a spatial composition approximately like this:

```cpp
enum class SpaceRole {
    PrimaryVoid,
    SecondaryVoid,
    Pocket,
    Throat,
    Shaft,
    Chasm,
    Overlook,
    Crossing,
    Reveal,
    Occluder,
    AnchorMass,
    DividerMass
};

struct SpatialNode {
    SpaceRole role;

    float scale;
    float openness;
    float verticality;
    float enclosure;
    float prominence;

    Vec3 preferredAxis;
    Range elevationRange;

    std::vector<Relation> relations;
};

enum class RelationType {
    Connects,
    Contains,
    Overlooks,
    Occludes,
    Wraps,
    Crosses,
    StacksAbove,
    VisibleFrom,
    RevealedBy
};

struct Relation {
    int otherNode;
    RelationType type;
    float strength;
};
```

The graph solver turns these relationships into positions and orientations. A later realization pass chooses appropriate SDF primitives and carving operators.

---

## The central shift

Your generator should not ask:

> What room do I put here?

It should ask:

> Does this part of the level need to enclose, expose, divide, frame, occlude, connect, overlook, wrap, compress, or reveal?

Those verbs can directly control density operations. They are abstract enough to work for a mine, alien ruin, natural cave, megastructure, or fantasy dungeon while still producing strongly authored spatial compositions.
