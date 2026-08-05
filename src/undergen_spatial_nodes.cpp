#include "undergen_spatial_nodes.h"

#include "density_grid.h"
#include <godot_cpp/classes/random_number_generator.hpp>
#include <godot_cpp/core/math.hpp>
#include <godot_cpp/variant/utility_functions.hpp>
#include <map>
#include <set>

namespace godot {

namespace {

Ref<UnderGenTopologyNodeData> make_topology_node(
        const String &p_id,
        UnderGenTopologyNodeData::SpaceRole p_role,
        int p_band,
        bool p_traversable,
        const PackedStringArray &p_tags = PackedStringArray()) {
    Ref<UnderGenTopologyNodeData> node;
    node.instantiate();
    node->set_id(p_id);
    node->set_role(p_role);
    node->set_elevation_band(p_band);
    node->set_traversable(p_traversable);
    node->set_gameplay_tags(p_tags);
    return node;
}

Ref<UnderGenTopologyEdgeData> make_topology_edge(
        const String &p_id,
        const String &p_from,
        const String &p_to,
        UnderGenTopologyEdgeData::RelationType p_relation,
        UnderGenTopologyEdgeData::TraversalType p_traversal,
        float p_width,
        float p_exposure,
        UnderGenTopologyEdgeData::DestinationVisibility p_visibility,
        bool p_required = true) {
    Ref<UnderGenTopologyEdgeData> edge;
    edge.instantiate();
    edge->set_id(p_id);
    edge->set_from_id(p_from);
    edge->set_to_id(p_to);
    edge->set_relation(p_relation);
    edge->set_traversal_type(p_traversal);
    edge->set_width_start(p_width);
    edge->set_width_end(p_width);
    edge->set_exposure(p_exposure);
    edge->set_enclosure(1.0f - p_exposure);
    edge->set_destination_visibility(p_visibility);
    edge->set_required(p_required);
    return edge;
}

UnderGenTopologyNodeData::SpaceRole infer_legacy_role(const String &p_symbol) {
    const String value = p_symbol.to_lower();
    if (value.contains("entry") || value.contains("start")) return UnderGenTopologyNodeData::ROLE_ENTRY;
    if (value.contains("exit")) return UnderGenTopologyNodeData::ROLE_EXIT;
    if (value.contains("boss")) return UnderGenTopologyNodeData::ROLE_BAY;
    if (value.contains("reveal")) return UnderGenTopologyNodeData::ROLE_REVEAL;
    if (value.contains("overlook")) return UnderGenTopologyNodeData::ROLE_OVERLOOK;
    if (value.contains("shaft")) return UnderGenTopologyNodeData::ROLE_SHAFT;
    if (value.contains("chasm")) return UnderGenTopologyNodeData::ROLE_CHASM;
    if (value.contains("pocket")) return UnderGenTopologyNodeData::ROLE_POCKET;
    return UnderGenTopologyNodeData::ROLE_SECONDARY_VOID;
}

String role_name(UnderGenTopologyNodeData::SpaceRole p_role) {
    switch (p_role) {
        case UnderGenTopologyNodeData::ROLE_ENTRY: return "entry";
        case UnderGenTopologyNodeData::ROLE_REVEAL: return "reveal";
        case UnderGenTopologyNodeData::ROLE_PRIMARY_VOID: return "primary_void";
        case UnderGenTopologyNodeData::ROLE_SECONDARY_VOID: return "secondary_void";
        case UnderGenTopologyNodeData::ROLE_POCKET: return "pocket";
        case UnderGenTopologyNodeData::ROLE_BAY: return "bay";
        case UnderGenTopologyNodeData::ROLE_THROAT: return "throat";
        case UnderGenTopologyNodeData::ROLE_SHAFT: return "shaft";
        case UnderGenTopologyNodeData::ROLE_CHASM: return "chasm";
        case UnderGenTopologyNodeData::ROLE_OVERLOOK: return "overlook";
        case UnderGenTopologyNodeData::ROLE_CROSSING: return "crossing";
        case UnderGenTopologyNodeData::ROLE_VISTA_TARGET: return "vista_target";
        case UnderGenTopologyNodeData::ROLE_ANCHOR_MASS: return "anchor_mass";
        case UnderGenTopologyNodeData::ROLE_OCCLUDER: return "occluder";
        case UnderGenTopologyNodeData::ROLE_DIVIDER_MASS: return "divider_mass";
        case UnderGenTopologyNodeData::ROLE_BOUNDARY_ROUTE: return "boundary_route";
        case UnderGenTopologyNodeData::ROLE_VERTICAL_CONNECTOR: return "vertical_connector";
        case UnderGenTopologyNodeData::ROLE_EXIT: return "exit";
        case UnderGenTopologyNodeData::ROLE_SLOT: return "slot";
        case UnderGenTopologyNodeData::ROLE_BOWL: return "bowl";
        case UnderGenTopologyNodeData::ROLE_DOME: return "dome";
        case UnderGenTopologyNodeData::ROLE_UNDERCROFT: return "undercroft";
        case UnderGenTopologyNodeData::ROLE_GALLERY: return "gallery";
        case UnderGenTopologyNodeData::ROLE_BUTTRESS: return "buttress";
        case UnderGenTopologyNodeData::ROLE_SPINE: return "spine";
        case UnderGenTopologyNodeData::ROLE_ISLAND: return "island";
        case UnderGenTopologyNodeData::ROLE_CANOPY: return "canopy";
        case UnderGenTopologyNodeData::ROLE_SCREEN: return "screen";
        case UnderGenTopologyNodeData::ROLE_DEAD_END: return "dead_end";
        case UnderGenTopologyNodeData::ROLE_LOOP_RETURN: return "loop_return";
        default: return "generic";
    }
}

bool is_mass_role(UnderGenTopologyNodeData::SpaceRole p_role) {
    return p_role == UnderGenTopologyNodeData::ROLE_ANCHOR_MASS ||
           p_role == UnderGenTopologyNodeData::ROLE_OCCLUDER ||
           p_role == UnderGenTopologyNodeData::ROLE_DIVIDER_MASS ||
           p_role == UnderGenTopologyNodeData::ROLE_BUTTRESS ||
           p_role == UnderGenTopologyNodeData::ROLE_SPINE ||
           p_role == UnderGenTopologyNodeData::ROLE_ISLAND ||
           p_role == UnderGenTopologyNodeData::ROLE_CANOPY ||
           p_role == UnderGenTopologyNodeData::ROLE_SCREEN;
}

float get_float(const Dictionary &p_dict, const StringName &p_key, float p_default) {
    return (float)p_dict.get(p_key, p_default);
}

Vector3 get_vec3(const Dictionary &p_dict, const StringName &p_key, const Vector3 &p_default) {
    return p_dict.get(p_key, p_default);
}

Vector3 clamp_inside(const Vector3 &p_position, const AABB &p_parent, const Vector3 &p_half_size, float p_margin) {
    Vector3 minp = p_parent.position + p_half_size + Vector3(p_margin, p_margin, p_margin);
    Vector3 maxp = p_parent.position + p_parent.size - p_half_size - Vector3(p_margin, p_margin, p_margin);
    return Vector3(
        Math::clamp(p_position.x, minp.x, maxp.x),
        Math::clamp(p_position.y, minp.y, maxp.y),
        Math::clamp(p_position.z, minp.z, maxp.z));
}

Ref<UnderGenPathState> make_path_state(float p_t, UnderGenPathState::SpatialState p_state,
        float p_width, float p_height, float p_exposure) {
    Ref<UnderGenPathState> state;
    state.instantiate();
    state->set_t(p_t);
    state->set_state(p_state);
    state->set_width(p_width);
    state->set_height(p_height);
    state->set_exposure(p_exposure);
    state->set_floor_flatness(p_state == UnderGenPathState::STATE_CREVICE ? 0.4f : 0.9f);
    if (p_state == UnderGenPathState::STATE_BOUNDARY || p_state == UnderGenPathState::STATE_LEDGE || p_state == UnderGenPathState::STATE_BALCONY) {
        state->set_right_wall(false);
    }
    return state;
}

Ref<UnderGenGeometryOperation> make_space_operation(const Ref<UnderGenEmbeddedSpace> &p_space,
        UnderGenGeometryOperation::OperationType p_type) {
    Ref<UnderGenGeometryOperation> op;
    op.instantiate();
    op->set_id(p_space->get_id() + "_shape");
    op->set_source_id(p_space->get_id());
    op->set_operation_type(p_type);
    op->set_transform(p_space->get_transform());
    op->set_size(p_space->get_size());
    op->set_zone_name(p_space->get_gameplay_tags().is_empty() ? role_name(p_space->get_role()) : p_space->get_gameplay_tags()[0]);
    switch (p_space->get_shape()) {
        case UnderGenEmbeddedSpace::SHAPE_BOX: op->set_primitive_type(UnderGenGeometryOperation::PRIMITIVE_BOX); break;
        case UnderGenEmbeddedSpace::SHAPE_CAPSULE: op->set_primitive_type(UnderGenGeometryOperation::PRIMITIVE_CAPSULE); break;
        default: op->set_primitive_type(UnderGenGeometryOperation::PRIMITIVE_ELLIPSOID); break;
    }
    return op;
}

Ref<UnderGenGeometryOperation> make_path_operation(const Ref<UnderGenEmbeddedPath> &p_path,
        UnderGenGeometryOperation::OperationType p_type, const PackedVector3Array &p_points,
        float p_radius, float p_height) {
    Ref<UnderGenGeometryOperation> op;
    op.instantiate();
    op->set_id(p_path->get_id() + "_" + String::num_int64((int)p_type));
    op->set_source_id(p_path->get_id());
    op->set_operation_type(p_type);
    op->set_primitive_type(UnderGenGeometryOperation::PRIMITIVE_SWEEP);
    op->set_points(p_points);
    op->set_radius(p_radius);
    op->set_height(p_height);
    Dictionary params = p_path->get_parameters();
    op->set_zone_name(params.get("zone_name", "route"));
    op->set_parameters(params);
    return op;
}

Vector3 sample_polyline(const PackedVector3Array &p_points, float p_t) {
    if (p_points.is_empty()) return Vector3();
    if (p_points.size() == 1 || p_t <= 0.0f) return p_points[0];
    if (p_t >= 1.0f) return p_points[p_points.size() - 1];
    float total = 0.0f;
    for (int i = 0; i < p_points.size() - 1; ++i) total += p_points[i].distance_to(p_points[i + 1]);
    if (total <= 0.001f) return p_points[0];
    float target = total * p_t;
    float walked = 0.0f;
    for (int i = 0; i < p_points.size() - 1; ++i) {
        float segment = p_points[i].distance_to(p_points[i + 1]);
        if (walked + segment >= target) return p_points[i].lerp(p_points[i + 1], (target - walked) / Math::max(segment, 0.001f));
        walked += segment;
    }
    return p_points[p_points.size() - 1];
}

} // namespace

void UnderGenTopologyBuilderNode::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_preset", "value"), &UnderGenTopologyBuilderNode::set_preset);
    ClassDB::bind_method(D_METHOD("get_preset"), &UnderGenTopologyBuilderNode::get_preset);
    ClassDB::bind_method(D_METHOD("set_prefer_legacy_input", "value"), &UnderGenTopologyBuilderNode::set_prefer_legacy_input);
    ClassDB::bind_method(D_METHOD("get_prefer_legacy_input"), &UnderGenTopologyBuilderNode::get_prefer_legacy_input);
    ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "preset", PROPERTY_HINT_RESOURCE_TYPE, "UnderGenGenerationPreset"), "set_preset", "get_preset");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "prefer_legacy_input"), "set_prefer_legacy_input", "get_prefer_legacy_input");
}

void UnderGenTopologyBuilderNode::set_preset(const Ref<UnderGenGenerationPreset> &v) { preset = v; emit_changed(); }
Ref<UnderGenGenerationPreset> UnderGenTopologyBuilderNode::get_preset() const { return preset; }
void UnderGenTopologyBuilderNode::set_prefer_legacy_input(bool v) { prefer_legacy_input = v; emit_changed(); }
bool UnderGenTopologyBuilderNode::get_prefer_legacy_input() const { return prefer_legacy_input; }

Dictionary UnderGenTopologyBuilderNode::get_pipeline_input_defaults(const Dictionary &global_inputs) const {
    Dictionary defaults;
    if (global_inputs.has(0)) defaults[0] = global_inputs[0];
    return defaults;
}

Ref<UnderGenSemanticGraph> UnderGenTopologyBuilderNode::_build_preset_graph(int64_t p_seed) const {
    Ref<UnderGenGenerationPreset> active = preset;
    if (active.is_null()) {
        active.instantiate();
        active->configure_builtin(UnderGenGenerationPreset::PRESET_LAYERED_CHASM);
    }

    Ref<UnderGenSemanticGraph> graph;
    graph.instantiate();
    graph->set_seed(p_seed);

    Dictionary metadata;
    metadata["preset_name"] = active->get_preset_name();
    metadata["world_size"] = active->get_world_size();
    metadata["primary_void_size"] = active->get_primary_void_size();
    metadata["anchor_scale"] = active->get_anchor_scale();
    metadata["anchor_eccentricity"] = active->get_anchor_eccentricity();
    metadata["elevation_band_count"] = active->get_elevation_band_count();
    metadata["elevation_band_spacing"] = active->get_elevation_band_spacing();
    metadata["route_width"] = active->get_route_width();
    metadata["route_height"] = active->get_route_height();
    metadata["openness"] = active->get_openness();
    metadata["verticality"] = active->get_verticality();
    metadata["exposure"] = active->get_exposure();
    metadata["noise_amplitude"] = active->get_noise_amplitude();
    graph->set_metadata(metadata);

    PackedStringArray entry_tags; entry_tags.append("entry"); entry_tags.append("player_start");
    PackedStringArray encounter_tags; encounter_tags.append("encounter");
    PackedStringArray exit_tags; exit_tags.append("exit"); exit_tags.append("portal");

    Ref<UnderGenTopologyNodeData> primary = make_topology_node("primary_void", UnderGenTopologyNodeData::ROLE_PRIMARY_VOID, 0, false);
    primary->set_openness(active->get_openness()); primary->set_verticality(active->get_verticality()); primary->set_prominence(1.0f);
    graph->add_node(primary);
    Ref<UnderGenTopologyNodeData> anchor = make_topology_node("anchor_mass", UnderGenTopologyNodeData::ROLE_ANCHOR_MASS, 0, false);
    anchor->set_prominence(1.0f); graph->add_node(anchor);
    Ref<UnderGenTopologyNodeData> lower_chasm = make_topology_node("lower_chasm", UnderGenTopologyNodeData::ROLE_CHASM, -1, false);
    lower_chasm->set_verticality(1.0f); lower_chasm->set_openness(active->get_openness()); lower_chasm->set_scale(1.15f); graph->add_node(lower_chasm);
    Ref<UnderGenTopologyNodeData> ceiling_pocket = make_topology_node("ceiling_pocket", UnderGenTopologyNodeData::ROLE_DOME, 1, false);
    ceiling_pocket->set_scale(0.8f); graph->add_node(ceiling_pocket);
    Ref<UnderGenTopologyNodeData> buttress = make_topology_node("secondary_buttress", UnderGenTopologyNodeData::ROLE_BUTTRESS, 0, false);
    buttress->set_prominence(0.65f); buttress->set_scale(0.7f); graph->add_node(buttress);
    graph->add_node(make_topology_node("entry_throat", UnderGenTopologyNodeData::ROLE_THROAT, 0, true, entry_tags));
    graph->add_node(make_topology_node("reveal", UnderGenTopologyNodeData::ROLE_REVEAL, 0, true));
    graph->add_node(make_topology_node("middle_band", UnderGenTopologyNodeData::ROLE_BOUNDARY_ROUTE, 0, true, encounter_tags));
    const bool has_upper_band = active->get_elevation_band_count() >= 2;
    const bool has_lower_band = active->get_elevation_band_count() >= 3;
    if (has_upper_band) graph->add_node(make_topology_node("upper_overlook", UnderGenTopologyNodeData::ROLE_OVERLOOK, 1, true));
    if (has_lower_band) graph->add_node(make_topology_node("lower_loop", UnderGenTopologyNodeData::ROLE_LOOP_RETURN, -1, true, encounter_tags));
    for (int i = 0; i < active->get_crossing_count(); ++i) {
        String crossing_id = i == 0 ? "crossing" : "crossing_" + String::num_int64(i + 1);
        Ref<UnderGenTopologyNodeData> crossing = make_topology_node(crossing_id, UnderGenTopologyNodeData::ROLE_CROSSING, has_upper_band ? 1 : 0, true);
        Dictionary crossing_params; crossing_params["crossing_index"] = i; crossing->set_parameters(crossing_params); graph->add_node(crossing);
    }
    graph->add_node(make_topology_node("exit_bay", UnderGenTopologyNodeData::ROLE_EXIT, 0, true, exit_tags));

    const float width = active->get_route_width();
    const float exposure = active->get_exposure();
    graph->add_edge(make_topology_edge("entry_to_reveal", "entry_throat", "reveal", UnderGenTopologyEdgeData::RELATION_CONNECTS, UnderGenTopologyEdgeData::TRAVERSAL_INTERIOR, width * 0.7f, 0.05f, UnderGenTopologyEdgeData::VISIBILITY_HIDDEN));
    graph->add_edge(make_topology_edge("reveal_to_middle", "reveal", "middle_band", UnderGenTopologyEdgeData::RELATION_CONNECTS, UnderGenTopologyEdgeData::TRAVERSAL_BOUNDARY, width, exposure, UnderGenTopologyEdgeData::VISIBILITY_PARTIAL));
    if (has_upper_band) graph->add_edge(make_topology_edge("reveal_to_upper", "reveal", "upper_overlook", UnderGenTopologyEdgeData::RELATION_CONNECTS, UnderGenTopologyEdgeData::TRAVERSAL_VERTICAL, width * 0.8f, exposure * 0.6f, UnderGenTopologyEdgeData::VISIBILITY_PARTIAL));
    String crossing_source = has_upper_band ? "upper_overlook" : "middle_band";
    for (int i = 0; i < active->get_crossing_count(); ++i) {
        String crossing_id = i == 0 ? "crossing" : "crossing_" + String::num_int64(i + 1);
        graph->add_edge(make_topology_edge("to_" + crossing_id, crossing_source, crossing_id, UnderGenTopologyEdgeData::RELATION_CROSSES, UnderGenTopologyEdgeData::TRAVERSAL_CROSSING, width * 0.8f, 1.0f, UnderGenTopologyEdgeData::VISIBILITY_CONTINUOUS));
        crossing_source = crossing_id;
    }
    graph->add_edge(make_topology_edge("crossing_route_to_exit", crossing_source, "exit_bay", UnderGenTopologyEdgeData::RELATION_CONNECTS, active->get_crossing_count() > 0 ? UnderGenTopologyEdgeData::TRAVERSAL_BOUNDARY : UnderGenTopologyEdgeData::TRAVERSAL_INTERIOR, width, exposure, UnderGenTopologyEdgeData::VISIBILITY_PARTIAL));
    if (has_lower_band) {
        graph->add_edge(make_topology_edge("middle_to_lower", "middle_band", "lower_loop", UnderGenTopologyEdgeData::RELATION_CONNECTS, UnderGenTopologyEdgeData::TRAVERSAL_VERTICAL, width * 0.75f, exposure * 0.5f, UnderGenTopologyEdgeData::VISIBILITY_PARTIAL));
        graph->add_edge(make_topology_edge("lower_to_exit", "lower_loop", "exit_bay", UnderGenTopologyEdgeData::RELATION_CONNECTS, UnderGenTopologyEdgeData::TRAVERSAL_BOUNDARY, width, exposure, UnderGenTopologyEdgeData::VISIBILITY_PARTIAL));
    }

    graph->add_edge(make_topology_edge("primary_contains_anchor", "primary_void", "anchor_mass", UnderGenTopologyEdgeData::RELATION_CONTAINS, UnderGenTopologyEdgeData::TRAVERSAL_NONE, width, 0.0f, UnderGenTopologyEdgeData::VISIBILITY_PARTIAL));
    graph->add_edge(make_topology_edge("primary_contains_chasm", "primary_void", "lower_chasm", UnderGenTopologyEdgeData::RELATION_CONTAINS, UnderGenTopologyEdgeData::TRAVERSAL_NONE, width, 1.0f, UnderGenTopologyEdgeData::VISIBILITY_PARTIAL));
    graph->add_edge(make_topology_edge("primary_contains_ceiling_pocket", "primary_void", "ceiling_pocket", UnderGenTopologyEdgeData::RELATION_CONTAINS, UnderGenTopologyEdgeData::TRAVERSAL_NONE, width, 0.4f, UnderGenTopologyEdgeData::VISIBILITY_PARTIAL));
    graph->add_edge(make_topology_edge("primary_contains_buttress", "primary_void", "secondary_buttress", UnderGenTopologyEdgeData::RELATION_CONTAINS, UnderGenTopologyEdgeData::TRAVERSAL_NONE, width, 0.0f, UnderGenTopologyEdgeData::VISIBILITY_PARTIAL));
    graph->add_edge(make_topology_edge("anchor_occludes_reveal", "anchor_mass", "reveal", UnderGenTopologyEdgeData::RELATION_OCCLUDES, UnderGenTopologyEdgeData::TRAVERSAL_NONE, width, 0.0f, UnderGenTopologyEdgeData::VISIBILITY_HIDDEN));
    graph->add_edge(make_topology_edge("middle_wraps_anchor", "middle_band", "anchor_mass", UnderGenTopologyEdgeData::RELATION_WRAPS, UnderGenTopologyEdgeData::TRAVERSAL_NONE, width, exposure, UnderGenTopologyEdgeData::VISIBILITY_PARTIAL));
    if (has_upper_band) {
        graph->add_edge(make_topology_edge("upper_overlooks_void", "upper_overlook", "primary_void", UnderGenTopologyEdgeData::RELATION_OVERLOOKS, UnderGenTopologyEdgeData::TRAVERSAL_NONE, width, 1.0f, UnderGenTopologyEdgeData::VISIBILITY_CONTINUOUS));
        graph->add_edge(make_topology_edge("upper_above_middle", "upper_overlook", "middle_band", UnderGenTopologyEdgeData::RELATION_STACKS_ABOVE, UnderGenTopologyEdgeData::TRAVERSAL_NONE, width, 0.0f, UnderGenTopologyEdgeData::VISIBILITY_CONTINUOUS));
    }
    if (has_lower_band) {
        graph->add_edge(make_topology_edge("middle_above_lower", "middle_band", "lower_loop", UnderGenTopologyEdgeData::RELATION_STACKS_ABOVE, UnderGenTopologyEdgeData::TRAVERSAL_NONE, width, 0.0f, UnderGenTopologyEdgeData::VISIBILITY_CONTINUOUS));
        graph->add_edge(make_topology_edge("lower_visible_from_middle", "lower_loop", "middle_band", UnderGenTopologyEdgeData::RELATION_VISIBLE_FROM, UnderGenTopologyEdgeData::TRAVERSAL_NONE, width, 1.0f, UnderGenTopologyEdgeData::VISIBILITY_PARTIAL));
    }

    for (int i = 0; i < active->get_secondary_bay_count(); ++i) {
        String id = "secondary_bay_" + String::num_int64(i + 1);
        PackedStringArray tags; tags.append("optional");
        Ref<UnderGenTopologyNodeData> bay = make_topology_node(id, i % 2 == 0 ? UnderGenTopologyNodeData::ROLE_BAY : UnderGenTopologyNodeData::ROLE_POCKET, (i % 3) - 1, true, tags);
        bay->set_scale(0.75f + 0.12f * (i % 3)); graph->add_node(bay);
        String source = i % 2 == 0 || !has_lower_band ? "middle_band" : "lower_loop";
        graph->add_edge(make_topology_edge("optional_to_" + id, source, id, UnderGenTopologyEdgeData::RELATION_CONNECTS, UnderGenTopologyEdgeData::TRAVERSAL_BOUNDARY, width * 0.75f, exposure * 0.5f, UnderGenTopologyEdgeData::VISIBILITY_HIDDEN, false));
    }

    return graph;
}

Ref<UnderGenSemanticGraph> UnderGenTopologyBuilderNode::_convert_legacy_graph(const Dictionary &p_graph, int64_t p_seed) const {
    Ref<UnderGenSemanticGraph> graph = _build_preset_graph(p_seed);
    // Retain only the macro forms from the preset; traversal comes from the grammar.
    TypedArray<UnderGenTopologyNodeData> macro_nodes;
    for (int i = 0; i < graph->get_nodes().size(); ++i) {
        Ref<UnderGenTopologyNodeData> n = graph->get_nodes()[i];
        if (n.is_valid() && !n->get_traversable()) macro_nodes.append(n);
    }
    graph->set_nodes(macro_nodes);
    TypedArray<UnderGenTopologyEdgeData> no_edges; graph->set_edges(no_edges);

    Array nodes = p_graph.get("nodes", Array());
    for (int i = 0; i < nodes.size(); ++i) {
        Dictionary legacy = nodes[i];
        String id = legacy.get("id", "legacy_" + String::num_int64(i));
        String symbol = legacy.get("symbol", legacy.get("type", "room"));
        PackedStringArray tags; tags.append(symbol.to_lower());
        Ref<UnderGenTopologyNodeData> node = make_topology_node(id, infer_legacy_role(symbol), 0, true, tags);
        Dictionary params;
        params["legacy_symbol"] = symbol;
        params["min_size"] = legacy.get("min_size", Vector3i(6, 5, 6));
        params["max_size"] = legacy.get("max_size", Vector3i(12, 8, 12));
        params["vox_path"] = legacy.get("vox_path", "");
        params["constraints"] = legacy.get("constraints", Dictionary());
        params["exclude_from_smoothing"] = legacy.get("exclude_from_smoothing", false);
        params["exclude_from_warping"] = legacy.get("exclude_from_warping", false);
        node->set_parameters(params);
        graph->add_node(node);
        graph->add_edge(make_topology_edge("primary_contains_" + id, "primary_void", id, UnderGenTopologyEdgeData::RELATION_CONTAINS, UnderGenTopologyEdgeData::TRAVERSAL_NONE, 4.0f, 0.0f, UnderGenTopologyEdgeData::VISIBILITY_PARTIAL));
    }

    Array edges = p_graph.get("edges", Array());
    for (int i = 0; i < edges.size(); ++i) {
        Dictionary legacy = edges[i];
        String type = String(legacy.get("type", "corridor")).to_lower();
        UnderGenTopologyEdgeData::TraversalType traversal = UnderGenTopologyEdgeData::TRAVERSAL_INTERIOR;
        if (type.contains("bridge") || type.contains("cross")) traversal = UnderGenTopologyEdgeData::TRAVERSAL_CROSSING;
        else if (type.contains("vertical") || type.contains("shaft")) traversal = UnderGenTopologyEdgeData::TRAVERSAL_VERTICAL;
        else if (type.contains("ledge") || type.contains("boundary")) traversal = UnderGenTopologyEdgeData::TRAVERSAL_BOUNDARY;
        Ref<UnderGenTopologyEdgeData> edge = make_topology_edge(
            "legacy_edge_" + String::num_int64(i), legacy.get("from", ""), legacy.get("to", ""),
            UnderGenTopologyEdgeData::RELATION_CONNECTS, traversal, 4.5f,
            traversal == UnderGenTopologyEdgeData::TRAVERSAL_BOUNDARY ? 0.7f : 0.2f,
            UnderGenTopologyEdgeData::VISIBILITY_PARTIAL);
        Dictionary params; params["legacy_type"] = type; edge->set_parameters(params); graph->add_edge(edge);
    }
    Dictionary metadata = graph->get_metadata(); metadata["legacy_source"] = true; graph->set_metadata(metadata);
    return graph;
}

void UnderGenTopologyBuilderNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    int64_t seed = inputs.get(0, (int64_t)12345);
    Variant legacy_variant = inputs.get(1, Variant());
    Ref<UnderGenSemanticGraph> graph;
    if (prefer_legacy_input && legacy_variant.get_type() == Variant::DICTIONARY && !((Dictionary)legacy_variant).is_empty()) graph = _convert_legacy_graph(legacy_variant, seed);
    else graph = _build_preset_graph(seed);
    PackedStringArray errors = graph->validate_graph();
    if (!errors.is_empty()) UtilityFunctions::printerr("UnderGenTopologyBuilderNode: generated topology has ", errors.size(), " validation error(s): ", errors);
    outputs[0] = graph;
}

void UnderGenSpatialEmbedderNode::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_solver_iterations", "value"), &UnderGenSpatialEmbedderNode::set_solver_iterations); ClassDB::bind_method(D_METHOD("get_solver_iterations"), &UnderGenSpatialEmbedderNode::get_solver_iterations);
    ClassDB::bind_method(D_METHOD("set_constraint_tolerance", "value"), &UnderGenSpatialEmbedderNode::set_constraint_tolerance); ClassDB::bind_method(D_METHOD("get_constraint_tolerance"), &UnderGenSpatialEmbedderNode::get_constraint_tolerance);
    ClassDB::bind_method(D_METHOD("set_preserve_vertical_bands", "value"), &UnderGenSpatialEmbedderNode::set_preserve_vertical_bands); ClassDB::bind_method(D_METHOD("get_preserve_vertical_bands"), &UnderGenSpatialEmbedderNode::get_preserve_vertical_bands);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "solver_iterations", PROPERTY_HINT_RANGE, "1,256,1"), "set_solver_iterations", "get_solver_iterations");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "constraint_tolerance", PROPERTY_HINT_RANGE, "0.001,1,0.001"), "set_constraint_tolerance", "get_constraint_tolerance");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "preserve_vertical_bands"), "set_preserve_vertical_bands", "get_preserve_vertical_bands");
}
void UnderGenSpatialEmbedderNode::set_solver_iterations(int v) { solver_iterations = Math::max(v, 1); emit_changed(); } int UnderGenSpatialEmbedderNode::get_solver_iterations() const { return solver_iterations; }
void UnderGenSpatialEmbedderNode::set_constraint_tolerance(float v) { constraint_tolerance = Math::max(v, 0.0001f); emit_changed(); } float UnderGenSpatialEmbedderNode::get_constraint_tolerance() const { return constraint_tolerance; }
void UnderGenSpatialEmbedderNode::set_preserve_vertical_bands(bool v) { preserve_vertical_bands = v; emit_changed(); } bool UnderGenSpatialEmbedderNode::get_preserve_vertical_bands() const { return preserve_vertical_bands; }
Dictionary UnderGenSpatialEmbedderNode::get_pipeline_input_defaults(const Dictionary &global_inputs) const { Dictionary result; if (global_inputs.has(0)) result[0] = global_inputs[0]; return result; }

Ref<UnderGenEmbeddedLayout> UnderGenSpatialEmbedderNode::_embed(const Ref<UnderGenSemanticGraph> &p_graph, int64_t p_seed) const {
    Ref<UnderGenEmbeddedLayout> layout; layout.instantiate(); layout->set_source_graph(p_graph);
    Dictionary metadata = p_graph->get_metadata();
    Vector3i world_i = metadata.get("world_size", Vector3i(192, 96, 192));
    Vector3 world = Vector3(world_i); Vector3 center = world * 0.5f;
    Vector3 primary_size = get_vec3(metadata, "primary_void_size", Vector3(118, 68, 92));
    float spacing = get_float(metadata, "elevation_band_spacing", 16.0f);
    float route_width = get_float(metadata, "route_width", 5.0f);
    float route_height = get_float(metadata, "route_height", 6.0f);
    float anchor_scale = get_float(metadata, "anchor_scale", 0.34f);
    float anchor_eccentricity = get_float(metadata, "anchor_eccentricity", 0.22f);

    Ref<RandomNumberGenerator> rng; rng.instantiate(); rng->set_seed(p_seed);
    int generic_index = 0; int bay_index = 0;
    TypedArray<UnderGenTopologyNodeData> nodes = p_graph->get_nodes();
    for (int i = 0; i < nodes.size(); ++i) {
        Ref<UnderGenTopologyNodeData> node = nodes[i]; if (node.is_null()) continue;
        Ref<UnderGenEmbeddedSpace> space; space.instantiate();
        space->set_id(node->get_id()); space->set_role(node->get_role()); space->set_gameplay_tags(node->get_gameplay_tags()); space->set_traversable(node->get_traversable());
        Dictionary space_params = node->get_parameters(); space_params["elevation_band"] = node->get_elevation_band(); space_params["elevation_range"] = node->get_elevation_range(); space->set_parameters(space_params);
        Vector3 position = center; Vector3 size(route_width * 2.0f, route_height * 1.5f, route_width * 2.0f);
        const float band_y = center.y + (preserve_vertical_bands ? node->get_elevation_band() * spacing : 0.0f);
        position.y = band_y;
        switch (node->get_role()) {
            case UnderGenTopologyNodeData::ROLE_PRIMARY_VOID:
                size = primary_size; position = center; space->set_shape(UnderGenEmbeddedSpace::SHAPE_ELLIPSOID); break;
            case UnderGenTopologyNodeData::ROLE_ANCHOR_MASS:
            case UnderGenTopologyNodeData::ROLE_OCCLUDER:
            case UnderGenTopologyNodeData::ROLE_DIVIDER_MASS:
            case UnderGenTopologyNodeData::ROLE_BUTTRESS:
            case UnderGenTopologyNodeData::ROLE_SPINE:
            case UnderGenTopologyNodeData::ROLE_ISLAND:
            case UnderGenTopologyNodeData::ROLE_CANOPY:
            case UnderGenTopologyNodeData::ROLE_SCREEN:
                size = Vector3(primary_size.x * anchor_scale, primary_size.y * 0.88f, primary_size.z * anchor_scale);
                position = center + Vector3(primary_size.x * anchor_eccentricity, 0, primary_size.z * 0.04f);
                if (node->get_role() == UnderGenTopologyNodeData::ROLE_BUTTRESS) { size *= Vector3(0.42f, 0.78f, 0.62f); position = center + Vector3(-primary_size.x * 0.34f, -primary_size.y * 0.06f, primary_size.z * 0.19f); }
                else if (node->get_role() == UnderGenTopologyNodeData::ROLE_CANOPY) { size *= Vector3(1.4f, 0.24f, 1.2f); position.y = center.y + primary_size.y * 0.34f; }
                else if (node->get_role() == UnderGenTopologyNodeData::ROLE_SCREEN) size *= Vector3(0.18f, 0.75f, 1.5f);
                space->set_shape(UnderGenEmbeddedSpace::SHAPE_MASS); break;
            case UnderGenTopologyNodeData::ROLE_CHASM:
                size = Vector3(primary_size.x * 0.56f, primary_size.y * 0.54f, primary_size.z * 0.58f) * node->get_scale();
                position = center + Vector3(primary_size.x * 0.08f, -primary_size.y * 0.30f, primary_size.z * 0.05f);
                space->set_shape(UnderGenEmbeddedSpace::SHAPE_ELLIPSOID); break;
            case UnderGenTopologyNodeData::ROLE_DOME:
                size = Vector3(primary_size.x * 0.42f, primary_size.y * 0.30f, primary_size.z * 0.38f) * node->get_scale();
                position = center + Vector3(-primary_size.x * 0.12f, primary_size.y * 0.38f, -primary_size.z * 0.12f);
                space->set_shape(UnderGenEmbeddedSpace::SHAPE_ELLIPSOID); break;
            case UnderGenTopologyNodeData::ROLE_ENTRY:
            case UnderGenTopologyNodeData::ROLE_THROAT:
                position += Vector3(-primary_size.x * 0.51f, 0, -primary_size.z * 0.22f);
                size = Vector3(route_width * 2.4f, route_height * 1.3f, route_width * 3.2f); space->set_shape(UnderGenEmbeddedSpace::SHAPE_CAPSULE); break;
            case UnderGenTopologyNodeData::ROLE_REVEAL:
                position += Vector3(-primary_size.x * 0.37f, 0, -primary_size.z * 0.17f);
                size = Vector3(route_width * 3.0f, route_height * 2.0f, route_width * 3.0f); break;
            case UnderGenTopologyNodeData::ROLE_OVERLOOK:
                position += Vector3(-primary_size.x * 0.25f, 0, -primary_size.z * 0.34f);
                size = Vector3(route_width * 3.0f, route_height * 1.5f, route_width * 2.2f); break;
            case UnderGenTopologyNodeData::ROLE_CROSSING:
                position += Vector3(primary_size.x * (0.27f - 0.10f * (int)node->get_parameters().get("crossing_index", 0)), 0, primary_size.z * (-0.08f + 0.18f * (int)node->get_parameters().get("crossing_index", 0)));
                size = Vector3(route_width * 2.0f, route_height * 1.2f, route_width * 2.0f); break;
            case UnderGenTopologyNodeData::ROLE_BOUNDARY_ROUTE:
            case UnderGenTopologyNodeData::ROLE_LOOP_RETURN:
                position += node->get_elevation_band() < 0
                    ? Vector3(-primary_size.x * 0.12f, 0, primary_size.z * 0.33f)
                    : Vector3(-primary_size.x * 0.30f, 0, primary_size.z * 0.22f);
                size = Vector3(route_width * 2.5f, route_height * 1.3f, route_width * 3.0f); break;
            case UnderGenTopologyNodeData::ROLE_EXIT:
                position += Vector3(primary_size.x * 0.46f, 0, primary_size.z * 0.20f);
                size = Vector3(route_width * 3.5f, route_height * 2.2f, route_width * 3.5f); break;
            case UnderGenTopologyNodeData::ROLE_BAY:
            case UnderGenTopologyNodeData::ROLE_POCKET: {
                float side = bay_index % 2 == 0 ? 1.0f : -1.0f;
                float z = ((bay_index % 3) - 1) * primary_size.z * 0.28f;
                position += Vector3(side * primary_size.x * 0.43f, 0, z);
                size = Vector3(route_width * (3.0f + node->get_scale()), route_height * (1.8f + node->get_scale() * 0.3f), route_width * (3.0f + node->get_scale()));
                bay_index++; break;
            }
            default: {
                float angle = (float)generic_index * 2.399963f;
                position += Vector3(Math::cos(angle) * primary_size.x * 0.32f, 0, Math::sin(angle) * primary_size.z * 0.32f);
                generic_index++;
                Dictionary params = node->get_parameters();
                if (params.has("min_size") && params.has("max_size")) {
                    Vector3 min_size = Vector3((Vector3i)params["min_size"]); Vector3 max_size = Vector3((Vector3i)params["max_size"]); size = (min_size + max_size) * 0.5f;
                }
                break;
            }
        }
        if (node->get_role() != UnderGenTopologyNodeData::ROLE_PRIMARY_VOID && !is_mass_role(node->get_role())) {
            position += Vector3(rng->randf_range(-1.5f, 1.5f), rng->randf_range(-0.5f, 0.5f), rng->randf_range(-1.5f, 1.5f));
        }
        space->set_position(position); space->set_size(size); layout->add_space(space);
    }

    Ref<UnderGenEmbeddedSpace> primary = layout->find_space("primary_void");
    Ref<UnderGenEmbeddedSpace> anchor = layout->find_space("anchor_mass");
    int satisfied = 0; PackedStringArray violations;
    TypedArray<UnderGenTopologyEdgeData> edges = p_graph->get_edges();
    for (int iteration = 0; iteration < solver_iterations; ++iteration) {
        bool changed = false;
        for (int i = 0; i < edges.size(); ++i) {
            Ref<UnderGenTopologyEdgeData> edge = edges[i]; if (edge.is_null() || edge->get_traversal_type() != UnderGenTopologyEdgeData::TRAVERSAL_NONE) continue;
            Ref<UnderGenEmbeddedSpace> from = layout->find_space(edge->get_from_id()); Ref<UnderGenEmbeddedSpace> to = layout->find_space(edge->get_to_id()); if (from.is_null() || to.is_null()) continue;
            if (edge->get_relation() == UnderGenTopologyEdgeData::RELATION_CONTAINS) {
                Vector3 clamped = clamp_inside(to->get_position(), from->get_bounds(), to->get_size() * 0.5f, 1.0f);
                if (clamped.distance_to(to->get_position()) > constraint_tolerance) { to->set_position(to->get_position().lerp(clamped, 0.65f)); changed = true; }
            } else if (edge->get_relation() == UnderGenTopologyEdgeData::RELATION_STACKS_ABOVE) {
                float required_y = to->get_position().y + spacing * 0.75f;
                if (from->get_position().y < required_y) { Vector3 pos = from->get_position(); pos.y = Math::lerp(pos.y, required_y, 0.65f); from->set_position(pos); changed = true; }
            } else if (edge->get_relation() == UnderGenTopologyEdgeData::RELATION_OVERLOOKS && from->get_position().y <= to->get_position().y) {
                Vector3 pos = from->get_position(); pos.y = to->get_position().y + spacing; from->set_position(pos); changed = true;
            }
        }
        if (!changed) break;
    }

    for (int i = 0; i < edges.size(); ++i) {
        Ref<UnderGenTopologyEdgeData> edge = edges[i]; if (edge.is_null()) continue;
        Ref<UnderGenEmbeddedSpace> from = layout->find_space(edge->get_from_id()); Ref<UnderGenEmbeddedSpace> to = layout->find_space(edge->get_to_id());
        if (from.is_null() || to.is_null()) { violations.append("Missing embedded endpoint for " + edge->get_id()); continue; }
        if (edge->get_traversal_type() == UnderGenTopologyEdgeData::TRAVERSAL_NONE) { satisfied++; continue; }

        Ref<UnderGenEmbeddedPath> path; path.instantiate(); path->set_id(edge->get_id()); path->set_from_id(edge->get_from_id()); path->set_to_id(edge->get_to_id()); path->set_traversal_type(edge->get_traversal_type());
        Vector3 a = from->get_position(), b = to->get_position(); PackedVector3Array points; points.append(a);
        Ref<UnderGenEmbeddedSpace> wrap_target;
        for (int relation_index = 0; relation_index < edges.size(); ++relation_index) {
            Ref<UnderGenTopologyEdgeData> relation = edges[relation_index]; if (relation.is_null() || relation->get_relation() != UnderGenTopologyEdgeData::RELATION_WRAPS) continue;
            if (relation->get_from_id() == edge->get_from_id() || relation->get_from_id() == edge->get_to_id()) { wrap_target = layout->find_space(relation->get_to_id()); break; }
        }
        if (edge->get_traversal_type() == UnderGenTopologyEdgeData::TRAVERSAL_VERTICAL) {
            points.append(Vector3(a.x, (a.y + b.y) * 0.5f, a.z)); points.append(Vector3(b.x, (a.y + b.y) * 0.5f, b.z));
        } else if (edge->get_traversal_type() == UnderGenTopologyEdgeData::TRAVERSAL_BOUNDARY && wrap_target.is_valid()) {
            Vector3 radial = a - wrap_target->get_position(); radial.y = 0; if (radial.length_squared() < 0.01f) radial = b - wrap_target->get_position(); radial.y = 0; if (radial.length_squared() < 0.01f) radial = Vector3(1, 0, 0); radial.normalize();
            float start_angle = Math::atan2(radial.z, radial.x); float wrap_radius = Math::max(wrap_target->get_size().x, wrap_target->get_size().z) * 0.62f + route_width;
            for (int arc = 1; arc <= 2; ++arc) { float angle = start_angle + Math::deg_to_rad(70.0f * arc); Vector3 control = wrap_target->get_position() + Vector3(Math::cos(angle) * wrap_radius, Math::lerp(a.y, b.y, (float)arc / 3.0f) - wrap_target->get_position().y, Math::sin(angle) * wrap_radius); points.append(control); }
        } else if (edge->get_traversal_type() == UnderGenTopologyEdgeData::TRAVERSAL_BOUNDARY && primary.is_valid()) {
            Vector3 midpoint = (a + b) * 0.5f; Vector3 radial = midpoint - primary->get_position(); radial.y = 0;
            if (radial.length_squared() < 0.01f) radial = Vector3(0, 0, 1); radial.normalize();
            Vector3 tangent(-radial.z, 0, radial.x);
            points.append(a.lerp(midpoint, 0.65f) + tangent * primary_size.x * 0.06f);
            points.append(midpoint + radial * primary_size.x * 0.10f);
            points.append(midpoint.lerp(b, 0.65f) - tangent * primary_size.x * 0.05f);
        } else if (edge->get_traversal_type() != UnderGenTopologyEdgeData::TRAVERSAL_CROSSING) {
            Vector3 midpoint = (a + b) * 0.5f; Vector3 direction = b - a; Vector3 bend(-direction.z, 0, direction.x);
            if (bend.length_squared() > 0.01f) bend = bend.normalized() * direction.length() * edge->get_curvature() * 0.18f;
            midpoint += bend;
            if (anchor.is_valid()) {
                Vector3 horizontal = midpoint - anchor->get_position(); horizontal.y = 0;
                float clearance = anchor->get_size().x * 0.65f + route_width;
                if (horizontal.length() < clearance) { if (horizontal.length_squared() < 0.01f) horizontal = Vector3(1, 0, 0); midpoint += horizontal.normalized() * (clearance - horizontal.length()); }
            }
            points.append(midpoint);
        }
        points.append(b); path->set_points(points);

        UnderGenPathState::SpatialState middle_state = UnderGenPathState::STATE_ENCLOSED;
        if (edge->get_traversal_type() == UnderGenTopologyEdgeData::TRAVERSAL_BOUNDARY) middle_state = edge->get_exposure() > 0.75f ? UnderGenPathState::STATE_BALCONY : UnderGenPathState::STATE_LEDGE;
        else if (edge->get_traversal_type() == UnderGenTopologyEdgeData::TRAVERSAL_CROSSING) middle_state = UnderGenPathState::STATE_BRIDGE;
        else if (edge->get_traversal_type() == UnderGenTopologyEdgeData::TRAVERSAL_VERTICAL) middle_state = UnderGenPathState::STATE_CREVICE;
        TypedArray<UnderGenPathState> states;
        states.append(make_path_state(0.0f, edge->get_destination_visibility() == UnderGenTopologyEdgeData::VISIBILITY_HIDDEN ? UnderGenPathState::STATE_COMPRESSED : UnderGenPathState::STATE_ENCLOSED, edge->get_width_start(), route_height, edge->get_exposure() * 0.25f));
        states.append(make_path_state(0.5f, middle_state, (edge->get_width_start() + edge->get_width_end()) * 0.5f, route_height, edge->get_exposure()));
        states.append(make_path_state(1.0f, UnderGenPathState::STATE_ENCLOSED, edge->get_width_end(), route_height, edge->get_exposure() * 0.25f));
        path->set_states(states);
        Dictionary params = edge->get_parameters(); params["radius"] = (edge->get_width_start() + edge->get_width_end()) * 0.25f; params["height"] = route_height; params["exposure"] = edge->get_exposure(); params["enclosure"] = edge->get_enclosure(); params["required"] = edge->get_required(); params["zone_name"] = from->get_gameplay_tags().is_empty() ? "route" : from->get_gameplay_tags()[0]; path->set_parameters(params);
        layout->add_path(path); satisfied++;
    }

    // Coarse semantic fields remain independent of density and can guide later geometry/detail passes.
    const Vector3i field_resolution(8, 4, 8);
    const AABB field_bounds(Vector3(), world);
    Ref<UnderGenSpatialField> field_refs[8];
    for (int type = 0; type < 8; ++type) { field_refs[type].instantiate(); field_refs[type]->initialize((UnderGenSpatialField::FieldType)type, field_bounds, field_resolution, 0.0f); }
    for (int z = 0; z < field_resolution.z; ++z) for (int y = 0; y < field_resolution.y; ++y) for (int x = 0; x < field_resolution.x; ++x) {
        Vector3 normalized((float)x / (field_resolution.x - 1), (float)y / (field_resolution.y - 1), (float)z / (field_resolution.z - 1));
        Vector3 sample = normalized * world;
        Vector3 q = (sample - center) / (primary_size * 0.5f); float void_distance = q.length();
        float openness = Math::clamp((1.25f - void_distance) * get_float(metadata, "openness", 0.7f), 0.0f, 1.0f);
        float exposure_field = Math::clamp(1.0f - Math::abs(void_distance - 0.82f) * 2.2f, 0.0f, 1.0f) * get_float(metadata, "exposure", 0.6f);
        float anchor_influence = 0.0f;
        if (anchor.is_valid()) { Vector3 aq = (sample - anchor->get_position()) / (anchor->get_size() * 0.5f); anchor_influence = Math::clamp(1.3f - aq.length(), 0.0f, 1.0f); }
        float route_influence = 0.0f;
        TypedArray<UnderGenEmbeddedPath> embedded_paths = layout->get_paths();
        for (int p = 0; p < embedded_paths.size(); ++p) { Ref<UnderGenEmbeddedPath> path = embedded_paths[p]; if (path.is_null()) continue; PackedVector3Array path_points = path->get_points(); for (int k = 0; k < path_points.size(); ++k) route_influence = Math::max(route_influence, Math::clamp(1.0f - sample.distance_to(path_points[k]) / 24.0f, 0.0f, 1.0f)); }
        Vector3i cell(x, y, z);
        field_refs[UnderGenSpatialField::FIELD_OPENNESS]->set_value(cell, openness);
        field_refs[UnderGenSpatialField::FIELD_VERTICALITY]->set_value(cell, Math::clamp(get_float(metadata, "verticality", 0.6f) * (0.75f + Math::abs(q.y) * 0.25f), 0.0f, 1.0f));
        field_refs[UnderGenSpatialField::FIELD_EXPOSURE]->set_value(cell, exposure_field);
        field_refs[UnderGenSpatialField::FIELD_ENCLOSURE]->set_value(cell, 1.0f - exposure_field);
        field_refs[UnderGenSpatialField::FIELD_OCCLUSION]->set_value(cell, anchor_influence);
        field_refs[UnderGenSpatialField::FIELD_PROMINENCE]->set_value(cell, Math::max(anchor_influence, exposure_field * 0.6f));
        field_refs[UnderGenSpatialField::FIELD_CONNECTIVITY_PRESSURE]->set_value(cell, route_influence);
        field_refs[UnderGenSpatialField::FIELD_SURFACE_SUITABILITY]->set_value(cell, exposure_field * Math::clamp(1.0f - Math::abs(q.y) * 0.45f, 0.0f, 1.0f));
    }
    for (int type = 0; type < 8; ++type) layout->add_field(field_refs[type]);

    Dictionary report; report["satisfied_constraints"] = satisfied; report["violations"] = violations; report["score"] = edges.is_empty() ? 1.0f : (float)satisfied / (float)edges.size(); report["solver_iterations"] = solver_iterations; report["seed"] = p_seed; report["semantic_field_count"] = 8; layout->set_constraint_report(report); layout->clear_dirty_regions();
    return layout;
}

void UnderGenSpatialEmbedderNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    int64_t seed = inputs.get(0, (int64_t)12345);
    Ref<UnderGenSemanticGraph> graph = inputs.get(1, Ref<UnderGenSemanticGraph>());
    if (graph.is_null()) { UtilityFunctions::printerr("UnderGenSpatialEmbedderNode: missing SemanticGraph on port 1"); return; }
    outputs[0] = _embed(graph, seed);
}

void UnderGenSpatialValidatorNode::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_minimum_positive_mass_ratio", "value"), &UnderGenSpatialValidatorNode::set_minimum_positive_mass_ratio); ClassDB::bind_method(D_METHOD("get_minimum_positive_mass_ratio"), &UnderGenSpatialValidatorNode::get_minimum_positive_mass_ratio);
    ClassDB::bind_method(D_METHOD("set_maximum_positive_mass_ratio", "value"), &UnderGenSpatialValidatorNode::set_maximum_positive_mass_ratio); ClassDB::bind_method(D_METHOD("get_maximum_positive_mass_ratio"), &UnderGenSpatialValidatorNode::get_maximum_positive_mass_ratio);
    ClassDB::bind_method(D_METHOD("set_maximum_crossing_ratio", "value"), &UnderGenSpatialValidatorNode::set_maximum_crossing_ratio); ClassDB::bind_method(D_METHOD("get_maximum_crossing_ratio"), &UnderGenSpatialValidatorNode::get_maximum_crossing_ratio);
    ClassDB::bind_method(D_METHOD("set_minimum_elevation_bands", "value"), &UnderGenSpatialValidatorNode::set_minimum_elevation_bands); ClassDB::bind_method(D_METHOD("get_minimum_elevation_bands"), &UnderGenSpatialValidatorNode::get_minimum_elevation_bands);
    ClassDB::bind_method(D_METHOD("set_reject_invalid_layout", "value"), &UnderGenSpatialValidatorNode::set_reject_invalid_layout); ClassDB::bind_method(D_METHOD("get_reject_invalid_layout"), &UnderGenSpatialValidatorNode::get_reject_invalid_layout);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "minimum_positive_mass_ratio", PROPERTY_HINT_RANGE, "0,1,0.01"), "set_minimum_positive_mass_ratio", "get_minimum_positive_mass_ratio");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "maximum_positive_mass_ratio", PROPERTY_HINT_RANGE, "0,1,0.01"), "set_maximum_positive_mass_ratio", "get_maximum_positive_mass_ratio");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "maximum_crossing_ratio", PROPERTY_HINT_RANGE, "0,1,0.01"), "set_maximum_crossing_ratio", "get_maximum_crossing_ratio");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "minimum_elevation_bands", PROPERTY_HINT_RANGE, "1,8,1"), "set_minimum_elevation_bands", "get_minimum_elevation_bands");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "reject_invalid_layout"), "set_reject_invalid_layout", "get_reject_invalid_layout");
}
void UnderGenSpatialValidatorNode::set_minimum_positive_mass_ratio(float v) { minimum_positive_mass_ratio = Math::clamp(v, 0.0f, 1.0f); emit_changed(); } float UnderGenSpatialValidatorNode::get_minimum_positive_mass_ratio() const { return minimum_positive_mass_ratio; }
void UnderGenSpatialValidatorNode::set_maximum_positive_mass_ratio(float v) { maximum_positive_mass_ratio = Math::clamp(v, 0.0f, 1.0f); emit_changed(); } float UnderGenSpatialValidatorNode::get_maximum_positive_mass_ratio() const { return maximum_positive_mass_ratio; }
void UnderGenSpatialValidatorNode::set_maximum_crossing_ratio(float v) { maximum_crossing_ratio = Math::clamp(v, 0.0f, 1.0f); emit_changed(); } float UnderGenSpatialValidatorNode::get_maximum_crossing_ratio() const { return maximum_crossing_ratio; }
void UnderGenSpatialValidatorNode::set_minimum_elevation_bands(int v) { minimum_elevation_bands = Math::clamp(v, 1, 8); emit_changed(); } int UnderGenSpatialValidatorNode::get_minimum_elevation_bands() const { return minimum_elevation_bands; }
void UnderGenSpatialValidatorNode::set_reject_invalid_layout(bool v) { reject_invalid_layout = v; emit_changed(); } bool UnderGenSpatialValidatorNode::get_reject_invalid_layout() const { return reject_invalid_layout; }

void UnderGenSpatialValidatorNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    Ref<UnderGenEmbeddedLayout> layout = inputs.get(0, Ref<UnderGenEmbeddedLayout>());
    if (layout.is_null()) { UtilityFunctions::printerr("UnderGenSpatialValidatorNode: missing EmbeddedLayout on port 0"); return; }
    Dictionary report = layout->validate_layout();
    PackedStringArray violations = report.get("violations", PackedStringArray());
    float primary_volume = 0.0f;
    float mass_volume = 0.0f;
    std::set<int> bands;
    TypedArray<UnderGenEmbeddedSpace> spaces = layout->get_spaces();
    for (int i = 0; i < spaces.size(); ++i) {
        Ref<UnderGenEmbeddedSpace> space = spaces[i]; if (space.is_null()) continue;
        Vector3 size = space->get_size(); float volume = size.x * size.y * size.z;
        const auto role = space->get_role();
        if (role == UnderGenTopologyNodeData::ROLE_PRIMARY_VOID) primary_volume += volume;
        if (is_mass_role(role)) mass_volume += volume;
        if (space->get_traversable()) bands.insert((int)space->get_parameters().get("elevation_band", 0));
    }
    int traversal_count = 0; int crossing_count = 0; int boundary_count = 0;
    TypedArray<UnderGenEmbeddedPath> paths = layout->get_paths();
    for (int i = 0; i < paths.size(); ++i) {
        Ref<UnderGenEmbeddedPath> path = paths[i]; if (path.is_null()) continue;
        traversal_count++;
        if (path->get_traversal_type() == UnderGenTopologyEdgeData::TRAVERSAL_CROSSING) crossing_count++;
        if (path->get_traversal_type() == UnderGenTopologyEdgeData::TRAVERSAL_BOUNDARY) boundary_count++;
    }
    const float mass_ratio = primary_volume > 0.001f ? mass_volume / primary_volume : 0.0f;
    const float crossing_ratio = traversal_count > 0 ? (float)crossing_count / (float)traversal_count : 0.0f;
    const float boundary_ratio = traversal_count > 0 ? (float)boundary_count / (float)traversal_count : 0.0f;
    if (mass_ratio < minimum_positive_mass_ratio || mass_ratio > maximum_positive_mass_ratio) violations.append("Positive-to-negative mass ratio is outside the configured range");
    if (crossing_ratio > maximum_crossing_ratio) violations.append("Cross-void traversal ratio exceeds the configured maximum");
    if ((int)bands.size() < minimum_elevation_bands) violations.append("Not enough distinct traversal elevation bands");
    report["positive_mass_ratio"] = mass_ratio;
    report["crossing_ratio"] = crossing_ratio;
    report["boundary_route_ratio"] = boundary_ratio;
    report["visible_elevation_band_count"] = (int)bands.size();
    report["violations"] = violations;
    report["valid"] = violations.is_empty();
    layout->set_constraint_report(report);
    outputs[1] = report;
    if (!reject_invalid_layout || violations.is_empty()) outputs[0] = layout;
}

void UnderGenGeometryPlannerNode::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_generate_boundary_ledges", "value"), &UnderGenGeometryPlannerNode::set_generate_boundary_ledges); ClassDB::bind_method(D_METHOD("get_generate_boundary_ledges"), &UnderGenGeometryPlannerNode::get_generate_boundary_ledges);
    ClassDB::bind_method(D_METHOD("set_generate_undercuts", "value"), &UnderGenGeometryPlannerNode::set_generate_undercuts); ClassDB::bind_method(D_METHOD("get_generate_undercuts"), &UnderGenGeometryPlannerNode::get_generate_undercuts);
    ClassDB::bind_method(D_METHOD("set_ledge_thickness", "value"), &UnderGenGeometryPlannerNode::set_ledge_thickness); ClassDB::bind_method(D_METHOD("get_ledge_thickness"), &UnderGenGeometryPlannerNode::get_ledge_thickness);
    ClassDB::bind_method(D_METHOD("build_plan", "layout"), &UnderGenGeometryPlannerNode::build_plan);
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "generate_boundary_ledges"), "set_generate_boundary_ledges", "get_generate_boundary_ledges"); ADD_PROPERTY(PropertyInfo(Variant::BOOL, "generate_undercuts"), "set_generate_undercuts", "get_generate_undercuts"); ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "ledge_thickness", PROPERTY_HINT_RANGE, "0.5,8,0.25"), "set_ledge_thickness", "get_ledge_thickness");
}
void UnderGenGeometryPlannerNode::set_generate_boundary_ledges(bool v) { generate_boundary_ledges = v; emit_changed(); } bool UnderGenGeometryPlannerNode::get_generate_boundary_ledges() const { return generate_boundary_ledges; }
void UnderGenGeometryPlannerNode::set_generate_undercuts(bool v) { generate_undercuts = v; emit_changed(); } bool UnderGenGeometryPlannerNode::get_generate_undercuts() const { return generate_undercuts; }
void UnderGenGeometryPlannerNode::set_ledge_thickness(float v) { ledge_thickness = Math::max(v, 0.25f); emit_changed(); } float UnderGenGeometryPlannerNode::get_ledge_thickness() const { return ledge_thickness; }

Ref<UnderGenGeometryPlan> UnderGenGeometryPlannerNode::_build_plan(const Ref<UnderGenEmbeddedLayout> &p_layout) const {
    Ref<UnderGenGeometryPlan> plan; plan.instantiate(); plan->set_source_layout(p_layout); plan->set_source_revision(p_layout->get_revision());
    Dictionary metadata; if (p_layout->get_source_graph().is_valid()) metadata = p_layout->get_source_graph()->get_metadata(); plan->set_metadata(metadata);
    TypedArray<UnderGenEmbeddedSpace> spaces = p_layout->get_spaces();
    // Ordered CSG: dominant void, preserved masses, secondary/local voids.
    for (int phase = 0; phase < 3; ++phase) {
        for (int i = 0; i < spaces.size(); ++i) {
            Ref<UnderGenEmbeddedSpace> space = spaces[i]; if (space.is_null()) continue;
            const auto role = space->get_role();
            bool primary = role == UnderGenTopologyNodeData::ROLE_PRIMARY_VOID;
            bool mass = is_mass_role(role);
            if ((phase == 0 && primary) || (phase == 1 && mass) || (phase == 2 && !primary && !mass)) {
                plan->add_operation(make_space_operation(space, mass ? UnderGenGeometryOperation::OP_ADD_MASS : UnderGenGeometryOperation::OP_SUBTRACT_VOID));
            }
        }
    }

    Ref<UnderGenSpatialField> openness_field = p_layout->find_field(UnderGenSpatialField::FIELD_OPENNESS);
    Ref<UnderGenSpatialField> exposure_field = p_layout->find_field(UnderGenSpatialField::FIELD_EXPOSURE);
    TypedArray<UnderGenEmbeddedPath> paths = p_layout->get_paths();
    for (int i = 0; i < paths.size(); ++i) {
        Ref<UnderGenEmbeddedPath> path = paths[i]; if (path.is_null()) continue;
        Dictionary params = path->get_parameters(); float radius = get_float(params, "radius", 2.5f); float height = get_float(params, "height", radius * 2.0f); float exposure = get_float(params, "exposure", 0.0f);
        PackedVector3Array points = path->get_points();
        if (points.size() < 2) continue;
        if (path->get_traversal_type() == UnderGenTopologyEdgeData::TRAVERSAL_CROSSING) {
            PackedVector3Array bridge_points = points; for (int p = 0; p < bridge_points.size(); ++p) bridge_points.set(p, bridge_points[p] - Vector3(0, ledge_thickness * 0.5f, 0));
            plan->add_operation(make_path_operation(path, UnderGenGeometryOperation::OP_BRIDGE, bridge_points, radius * 0.85f, ledge_thickness));
        } else if (path->get_traversal_type() == UnderGenTopologyEdgeData::TRAVERSAL_BOUNDARY && generate_boundary_ledges) {
            PackedVector3Array ledge_points = points; for (int p = 0; p < ledge_points.size(); ++p) ledge_points.set(p, ledge_points[p] - Vector3(0, ledge_thickness * 0.5f, 0));
            plan->add_operation(make_path_operation(path, UnderGenGeometryOperation::OP_LEDGE, ledge_points, radius * 1.15f, ledge_thickness));
            if (generate_undercuts && exposure > 0.45f) {
                PackedVector3Array undercut = points; for (int p = 0; p < undercut.size(); ++p) undercut.set(p, undercut[p] - Vector3(0, ledge_thickness + radius * 0.6f, 0));
                plan->add_operation(make_path_operation(path, UnderGenGeometryOperation::OP_UNDERCUT, undercut, radius * 0.75f, radius));
            }
        }
        TypedArray<UnderGenPathState> states = path->get_states();
        if (states.size() < 2) {
            plan->add_operation(make_path_operation(path, UnderGenGeometryOperation::OP_ROUTE_CLEARANCE, points, radius, height));
            continue;
        }
        for (int s = 0; s < states.size() - 1; ++s) {
            Ref<UnderGenPathState> from_state = states[s]; Ref<UnderGenPathState> to_state = states[s + 1];
            if (from_state.is_null() || to_state.is_null()) continue;
            PackedVector3Array segment;
            Vector3 start = sample_polyline(points, from_state->get_t());
            Vector3 middle = sample_polyline(points, (from_state->get_t() + to_state->get_t()) * 0.5f);
            Vector3 end = sample_polyline(points, to_state->get_t());
            const float segment_height = (from_state->get_height() + to_state->get_height()) * 0.5f;
            const float local_openness = openness_field.is_valid() ? openness_field->sample(middle) : 0.5f;
            const float local_exposure = exposure_field.is_valid() ? exposure_field->sample(middle) : (from_state->get_exposure() + to_state->get_exposure()) * 0.5f;
            const float vertical_offset = path->get_traversal_type() == UnderGenTopologyEdgeData::TRAVERSAL_BOUNDARY ? segment_height * 0.32f : segment_height * 0.18f;
            segment.append(start + Vector3(0, vertical_offset, 0)); segment.append(middle + Vector3(0, vertical_offset, 0)); segment.append(end + Vector3(0, vertical_offset, 0));
            Ref<UnderGenGeometryOperation> clearance = make_path_operation(path, UnderGenGeometryOperation::OP_ROUTE_CLEARANCE, segment, (from_state->get_width() + to_state->get_width()) * 0.25f * Math::lerp(0.88f, 1.18f, local_openness), segment_height);
            clearance->set_id(path->get_id() + "_state_" + String::num_int64(s));
            Dictionary state_params = clearance->get_parameters();
            state_params["state_from"] = (int)from_state->get_state(); state_params["state_to"] = (int)to_state->get_state();
            state_params["left_wall"] = from_state->get_left_wall() && to_state->get_left_wall(); state_params["right_wall"] = from_state->get_right_wall() && to_state->get_right_wall();
            state_params["floor_flatness"] = (from_state->get_floor_flatness() + to_state->get_floor_flatness()) * 0.5f;
            state_params["lateral_bias"] = (from_state->get_lateral_bias() + to_state->get_lateral_bias()) * 0.5f;
            state_params["local_noise_scale"] = (from_state->get_local_noise_scale() + to_state->get_local_noise_scale()) * 0.5f;
            state_params["sampled_openness"] = local_openness; state_params["sampled_exposure"] = local_exposure;
            clearance->set_parameters(state_params);
            plan->add_operation(clearance);
        }
    }
    return plan;
}

Ref<UnderGenGeometryPlan> UnderGenGeometryPlannerNode::build_plan(const Ref<UnderGenEmbeddedLayout> &p_layout) const {
    if (p_layout.is_null()) return Ref<UnderGenGeometryPlan>();
    return _build_plan(p_layout);
}

void UnderGenGeometryPlannerNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    Ref<UnderGenEmbeddedLayout> layout = inputs.get(0, Ref<UnderGenEmbeddedLayout>());
    if (layout.is_null()) { UtilityFunctions::printerr("UnderGenGeometryPlannerNode: missing EmbeddedLayout on port 0"); return; }
    outputs[0] = _build_plan(layout);
}

void UnderGenGeometryRealizerNode::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_surface_threshold", "value"), &UnderGenGeometryRealizerNode::set_surface_threshold); ClassDB::bind_method(D_METHOD("get_surface_threshold"), &UnderGenGeometryRealizerNode::get_surface_threshold);
    ClassDB::bind_method(D_METHOD("set_voxel_size", "value"), &UnderGenGeometryRealizerNode::set_voxel_size); ClassDB::bind_method(D_METHOD("get_voxel_size"), &UnderGenGeometryRealizerNode::get_voxel_size);
    ClassDB::bind_method(D_METHOD("set_retain_spatial_data", "value"), &UnderGenGeometryRealizerNode::set_retain_spatial_data); ClassDB::bind_method(D_METHOD("get_retain_spatial_data"), &UnderGenGeometryRealizerNode::get_retain_spatial_data);
    ClassDB::bind_method(D_METHOD("realize_plan", "plan"), &UnderGenGeometryRealizerNode::realize_plan);
    ClassDB::bind_method(D_METHOD("rebuild_dirty_regions", "context", "plan", "dirty_regions"), &UnderGenGeometryRealizerNode::rebuild_dirty_regions);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "surface_threshold", PROPERTY_HINT_RANGE, "-1,1,0.01"), "set_surface_threshold", "get_surface_threshold"); ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "voxel_size", PROPERTY_HINT_RANGE, "0.05,8,0.05"), "set_voxel_size", "get_voxel_size"); ADD_PROPERTY(PropertyInfo(Variant::BOOL, "retain_spatial_data"), "set_retain_spatial_data", "get_retain_spatial_data");
}
void UnderGenGeometryRealizerNode::set_surface_threshold(float v) { surface_threshold = v; emit_changed(); } float UnderGenGeometryRealizerNode::get_surface_threshold() const { return surface_threshold; }
void UnderGenGeometryRealizerNode::set_voxel_size(float v) { voxel_size = Math::max(v, 0.01f); emit_changed(); } float UnderGenGeometryRealizerNode::get_voxel_size() const { return voxel_size; }
void UnderGenGeometryRealizerNode::set_retain_spatial_data(bool v) { retain_spatial_data = v; emit_changed(); } bool UnderGenGeometryRealizerNode::get_retain_spatial_data() const { return retain_spatial_data; }

void UnderGenGeometryRealizerNode::_apply_primitive(const Ref<DensityGrid> &p_grid, const Ref<UnderGenGeometryOperation> &op, bool p_add_mass, bool p_use_clip, const AABB &p_clip) const {
    AABB bounds = op->get_bounds(); Vector3i dims = p_grid->get_grid_dimensions();
    if (p_use_clip && !bounds.intersects(p_clip)) return;
    if (p_use_clip) bounds = bounds.intersection(p_clip);
    Vector3i minp(Math::max(0, (int)Math::floor(bounds.position.x)), Math::max(0, (int)Math::floor(bounds.position.y)), Math::max(0, (int)Math::floor(bounds.position.z)));
    Vector3 end = bounds.position + bounds.size;
    Vector3i maxp(Math::min(dims.x - 1, (int)Math::ceil(end.x)), Math::min(dims.y - 1, (int)Math::ceil(end.y)), Math::min(dims.z - 1, (int)Math::ceil(end.z)));
    Transform3D inverse = op->get_transform().affine_inverse(); Vector3 half = op->get_size() * 0.5f; half.x = Math::max(half.x, 0.01f); half.y = Math::max(half.y, 0.01f); half.z = Math::max(half.z, 0.01f);
    int zone = !p_add_mass && !op->get_zone_name().is_empty() ? p_grid->register_zone_name(op->get_zone_name()) : 0;
    for (int x = minp.x; x <= maxp.x; ++x) for (int y = minp.y; y <= maxp.y; ++y) for (int z = minp.z; z <= maxp.z; ++z) {
        Vector3 local = inverse.xform(Vector3(x, y, z)); bool inside = false;
        if (op->get_primitive_type() == UnderGenGeometryOperation::PRIMITIVE_BOX) inside = Math::abs(local.x) <= half.x && Math::abs(local.y) <= half.y && Math::abs(local.z) <= half.z;
        else {
            Vector3 q(local.x / half.x, local.y / half.y, local.z / half.z); inside = q.length_squared() <= 1.0f;
        }
        if (inside) { Vector3i cell(x, y, z); p_grid->set_cell_fast(cell, p_add_mass ? 1.0f : -1.0f); if (zone > 0 && !p_add_mass) p_grid->set_zone_at_fast(cell, zone); if (op->get_material_id() > 0) p_grid->set_material_id_fast(cell, op->get_material_id()); }
    }
}

void UnderGenGeometryRealizerNode::_apply_sweep(const Ref<DensityGrid> &p_grid, const Ref<UnderGenGeometryOperation> &op, bool p_add_mass, bool p_use_clip, const AABB &p_clip) const {
    PackedVector3Array points = op->get_points(); if (points.is_empty()) return;
    Vector3i dims = p_grid->get_grid_dimensions(); float radius = op->get_radius(); float half_height = Math::max(op->get_height() * 0.5f, 0.25f);
    Dictionary op_params = op->get_parameters();
    const bool left_wall = op_params.get("left_wall", true); const bool right_wall = op_params.get("right_wall", true);
    const float lateral_bias = get_float(op_params, "lateral_bias", 0.0f);
    const float open_extension = (!left_wall || !right_wall) ? 1.35f : 1.0f;
    int zone = !p_add_mass && !op->get_zone_name().is_empty() ? p_grid->register_zone_name(op->get_zone_name()) : 0;
    auto stamp = [&](const Vector3 &raw_center, const Vector3 &raw_tangent) {
        Vector3 tangent(raw_tangent.x, 0, raw_tangent.z); if (tangent.length_squared() < 0.001f) tangent = Vector3(0, 0, 1); else tangent.normalize();
        Vector3 side(-tangent.z, 0, tangent.x);
        Vector3 center = raw_center + side * lateral_bias * radius;
        float stamp_radius = radius * open_extension;
        Vector3i minp(Math::max(0, (int)Math::floor(center.x - stamp_radius)), Math::max(0, (int)Math::floor(center.y - half_height)), Math::max(0, (int)Math::floor(center.z - stamp_radius)));
        Vector3i maxp(Math::min(dims.x - 1, (int)Math::ceil(center.x + stamp_radius)), Math::min(dims.y - 1, (int)Math::ceil(center.y + half_height)), Math::min(dims.z - 1, (int)Math::ceil(center.z + stamp_radius)));
        if (p_use_clip) {
            Vector3 clip_end = p_clip.position + p_clip.size;
            minp.x = Math::max(minp.x, (int)Math::floor(p_clip.position.x)); minp.y = Math::max(minp.y, (int)Math::floor(p_clip.position.y)); minp.z = Math::max(minp.z, (int)Math::floor(p_clip.position.z));
            maxp.x = Math::min(maxp.x, (int)Math::ceil(clip_end.x)); maxp.y = Math::min(maxp.y, (int)Math::ceil(clip_end.y)); maxp.z = Math::min(maxp.z, (int)Math::ceil(clip_end.z));
        }
        for (int x = minp.x; x <= maxp.x; ++x) for (int y = minp.y; y <= maxp.y; ++y) for (int z = minp.z; z <= maxp.z; ++z) {
            Vector3 delta = Vector3(x, y, z) - center; float side_distance = delta.dot(side); float lateral_radius = radius;
            if ((!right_wall && side_distance > 0.0f) || (!left_wall && side_distance < 0.0f)) lateral_radius *= open_extension;
            float normalized = (delta.x * delta.x + delta.z * delta.z) / (lateral_radius * lateral_radius) + (delta.y * delta.y) / (half_height * half_height);
            if (normalized <= 1.0f) { Vector3i cell(x, y, z); p_grid->set_cell_fast(cell, p_add_mass ? 1.0f : -1.0f); if (zone > 0 && !p_add_mass) p_grid->set_zone_at_fast(cell, zone); if (op->get_material_id() > 0) p_grid->set_material_id_fast(cell, op->get_material_id()); }
        }
    };
    Vector3 first_tangent = points.size() > 1 ? points[1] - points[0] : Vector3(0, 0, 1);
    stamp(points[0], first_tangent);
    for (int i = 0; i < points.size() - 1; ++i) {
        Vector3 a = points[i], b = points[i + 1]; float length = a.distance_to(b); int steps = Math::max(1, (int)Math::ceil(length * 1.25f));
        Vector3 tangent = b - a;
        for (int step = 1; step <= steps; ++step) stamp(a.lerp(b, (float)step / (float)steps), tangent);
    }
}

void UnderGenGeometryRealizerNode::_apply_operation(const Ref<DensityGrid> &p_grid, const Ref<UnderGenGeometryOperation> &op, bool p_use_clip, const AABB &p_clip) const {
    bool add_mass = op->get_operation_type() == UnderGenGeometryOperation::OP_ADD_MASS || op->get_operation_type() == UnderGenGeometryOperation::OP_LEDGE || op->get_operation_type() == UnderGenGeometryOperation::OP_BRIDGE || op->get_operation_type() == UnderGenGeometryOperation::OP_TERRACE;
    if (p_use_clip && !op->get_bounds().intersects(p_clip)) return;
    if (op->get_primitive_type() == UnderGenGeometryOperation::PRIMITIVE_SWEEP) _apply_sweep(p_grid, op, add_mass, p_use_clip, p_clip); else _apply_primitive(p_grid, op, add_mass, p_use_clip, p_clip);
}

Dictionary UnderGenGeometryRealizerNode::_update_context_from_plan(const Dictionary &p_context, const Ref<UnderGenGeometryPlan> &plan, const Ref<DensityGrid> &grid) const {
    Dictionary context = p_context.duplicate();
    context["grid"] = grid;
    context["seed"] = plan->get_source_layout().is_valid() && plan->get_source_layout()->get_source_graph().is_valid() ? plan->get_source_layout()->get_source_graph()->get_seed() : (int64_t)12345;
    Array rooms; Array edges;
    if (plan->get_source_layout().is_valid()) {
        Ref<UnderGenEmbeddedLayout> layout = plan->get_source_layout();
        layout->validate_layout();
        TypedArray<UnderGenEmbeddedSpace> spaces = layout->get_spaces();
        for (int i = 0; i < spaces.size(); ++i) {
            Ref<UnderGenEmbeddedSpace> space = spaces[i]; if (space.is_null() || !space->get_traversable()) continue;
            Dictionary room; Vector3 size = space->get_size(); Vector3 position = space->get_position() - size * 0.5f;
            room["id"] = space->get_id(); room["type"] = space->get_gameplay_tags().is_empty() ? role_name(space->get_role()) : space->get_gameplay_tags()[0]; room["position"] = Vector3i(position); room["size"] = Vector3i(size); room["center"] = space->get_position(); room["spatial_role"] = role_name(space->get_role()); room["gameplay_tags"] = space->get_gameplay_tags(); rooms.append(room);
            Dictionary parameters = space->get_parameters();
            room["spatial_parameters"] = parameters.duplicate();
            room["vox_path"] = parameters.get("vox_path", "");
            room["constraints"] = parameters.get("constraints", Dictionary());
            room["exclude_from_smoothing"] = parameters.get("exclude_from_smoothing", false);
            room["exclude_from_warping"] = parameters.get("exclude_from_warping", false);
            rooms[rooms.size() - 1] = room;
        }
        Ref<UnderGenSemanticGraph> graph = layout->get_source_graph();
        if (graph.is_valid()) {
            TypedArray<UnderGenTopologyEdgeData> topology_edges = graph->get_edges();
            for (int i = 0; i < topology_edges.size(); ++i) { Ref<UnderGenTopologyEdgeData> edge = topology_edges[i]; if (edge.is_null() || edge->get_traversal_type() == UnderGenTopologyEdgeData::TRAVERSAL_NONE) continue; Dictionary e; e["from"] = edge->get_from_id(); e["to"] = edge->get_to_id(); e["type"] = (int)edge->get_traversal_type(); e["exposure"] = edge->get_exposure(); e["enclosure"] = edge->get_enclosure(); edges.append(e); }
        }
        if (retain_spatial_data) { context["semantic_graph"] = graph; context["embedded_layout"] = layout; context["geometry_plan"] = plan; context["constraint_report"] = layout->get_constraint_report(); }
    }
    context["rooms"] = rooms; context["edges"] = edges;
    return context;
}

Dictionary UnderGenGeometryRealizerNode::realize_plan(const Ref<UnderGenGeometryPlan> &plan) const {
    if (plan.is_null()) return Dictionary();
    Dictionary metadata = plan->get_metadata(); Vector3i world_size = metadata.get("world_size", Vector3i(192, 96, 192));
    Ref<DensityGrid> grid; grid.instantiate(); grid->initialize_grid(world_size.x, world_size.y, world_size.z, 1.0f); grid->set_surface_threshold(surface_threshold);
    TypedArray<UnderGenGeometryOperation> operations = plan->get_operations();
    for (int i = 0; i < operations.size(); ++i) { Ref<UnderGenGeometryOperation> op = operations[i]; if (op.is_valid()) _apply_operation(grid, op); }
    return _update_context_from_plan(Dictionary(), plan, grid);
}

Dictionary UnderGenGeometryRealizerNode::rebuild_dirty_regions(const Dictionary &p_context, const Ref<UnderGenGeometryPlan> &plan, const Array &p_dirty_regions) const {
    if (plan.is_null() || p_dirty_regions.is_empty()) return p_context;
    Ref<DensityGrid> grid = p_context.get("grid", Ref<DensityGrid>());
    if (grid.is_null()) return realize_plan(plan);
    AABB dirty;
    bool has_dirty = false;
    for (int i = 0; i < p_dirty_regions.size(); ++i) {
        if (p_dirty_regions[i].get_type() != Variant::AABB) continue;
        AABB region = p_dirty_regions[i];
        if (region.size == Vector3()) continue;
        dirty = has_dirty ? dirty.merge(region) : region;
        has_dirty = true;
    }
    if (!has_dirty) return _update_context_from_plan(p_context, plan, grid);
    dirty = dirty.grow(2.0f);
    Vector3i dims = grid->get_grid_dimensions();
    Vector3 end = dirty.position + dirty.size;
    Vector3i minp(Math::max(0, (int)Math::floor(dirty.position.x)), Math::max(0, (int)Math::floor(dirty.position.y)), Math::max(0, (int)Math::floor(dirty.position.z)));
    Vector3i maxp(Math::min(dims.x - 1, (int)Math::ceil(end.x)), Math::min(dims.y - 1, (int)Math::ceil(end.y)), Math::min(dims.z - 1, (int)Math::ceil(end.z)));
    PackedFloat32Array &density_data = grid->get_density_data_rw(); PackedByteArray &material_data = grid->get_material_data_rw(); PackedInt32Array &zone_data = grid->get_zone_data_rw();
    float *density_ptr = density_data.ptrw(); uint8_t *material_ptr = material_data.ptrw(); int32_t *zone_ptr = zone_data.ptrw();
    for (int z = minp.z; z <= maxp.z; ++z) for (int y = minp.y; y <= maxp.y; ++y) for (int x = minp.x; x <= maxp.x; ++x) {
        int index = grid->get_index_unchecked(x, y, z); density_ptr[index] = 1.0f; material_ptr[index] = 0; zone_ptr[index] = 0;
    }
    grid->mark_flat_cache_dirty();
    grid->clear_hermite_edges();
    TypedArray<UnderGenGeometryOperation> operations = plan->get_operations();
    for (int i = 0; i < operations.size(); ++i) { Ref<UnderGenGeometryOperation> op = operations[i]; if (op.is_valid()) _apply_operation(grid, op, true, dirty); }
    return _update_context_from_plan(p_context, plan, grid);
}

void UnderGenGeometryRealizerNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    Ref<UnderGenGeometryPlan> plan = inputs.get(0, Ref<UnderGenGeometryPlan>());
    if (plan.is_null()) { UtilityFunctions::printerr("UnderGenGeometryRealizerNode: missing GeometryPlan on port 0"); return; }
    outputs[0] = realize_plan(plan);
    outputs[1] = plan;
}

void UnderGenGameplayMarkerNode::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_emit_player_start", "value"), &UnderGenGameplayMarkerNode::set_emit_player_start);
    ClassDB::bind_method(D_METHOD("get_emit_player_start"), &UnderGenGameplayMarkerNode::get_emit_player_start);
    ClassDB::bind_method(D_METHOD("set_emit_exit_portals", "value"), &UnderGenGameplayMarkerNode::set_emit_exit_portals);
    ClassDB::bind_method(D_METHOD("get_emit_exit_portals"), &UnderGenGameplayMarkerNode::get_emit_exit_portals);
    ClassDB::bind_method(D_METHOD("set_emit_encounters", "value"), &UnderGenGameplayMarkerNode::set_emit_encounters);
    ClassDB::bind_method(D_METHOD("get_emit_encounters"), &UnderGenGameplayMarkerNode::get_emit_encounters);
    ClassDB::bind_method(D_METHOD("set_replace_authored_route_markers", "value"), &UnderGenGameplayMarkerNode::set_replace_authored_route_markers);
    ClassDB::bind_method(D_METHOD("get_replace_authored_route_markers"), &UnderGenGameplayMarkerNode::get_replace_authored_route_markers);
    ClassDB::bind_method(D_METHOD("set_floor_offset", "value"), &UnderGenGameplayMarkerNode::set_floor_offset);
    ClassDB::bind_method(D_METHOD("get_floor_offset"), &UnderGenGameplayMarkerNode::get_floor_offset);
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "emit_player_start"), "set_emit_player_start", "get_emit_player_start");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "emit_exit_portals"), "set_emit_exit_portals", "get_emit_exit_portals");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "emit_encounters"), "set_emit_encounters", "get_emit_encounters");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "replace_authored_route_markers"), "set_replace_authored_route_markers", "get_replace_authored_route_markers");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "floor_offset", PROPERTY_HINT_RANGE, "0,8,0.1"), "set_floor_offset", "get_floor_offset");
}

void UnderGenGameplayMarkerNode::set_emit_player_start(bool v) { emit_player_start = v; emit_changed(); }
bool UnderGenGameplayMarkerNode::get_emit_player_start() const { return emit_player_start; }
void UnderGenGameplayMarkerNode::set_emit_exit_portals(bool v) { emit_exit_portals = v; emit_changed(); }
bool UnderGenGameplayMarkerNode::get_emit_exit_portals() const { return emit_exit_portals; }
void UnderGenGameplayMarkerNode::set_emit_encounters(bool v) { emit_encounters = v; emit_changed(); }
bool UnderGenGameplayMarkerNode::get_emit_encounters() const { return emit_encounters; }
void UnderGenGameplayMarkerNode::set_replace_authored_route_markers(bool v) { replace_authored_route_markers = v; emit_changed(); }
bool UnderGenGameplayMarkerNode::get_replace_authored_route_markers() const { return replace_authored_route_markers; }
void UnderGenGameplayMarkerNode::set_floor_offset(float v) { floor_offset = Math::max(v, 0.0f); emit_changed(); }
float UnderGenGameplayMarkerNode::get_floor_offset() const { return floor_offset; }

void UnderGenGameplayMarkerNode::_execute(const Dictionary &inputs, Dictionary &outputs) {
    Dictionary context = inputs.get(0, Dictionary());
    if (context.is_empty()) return;

    Array preserved_spawns;
    Array existing_spawns = context.get("vox_spawns", Array());
    for (int i = 0; i < existing_spawns.size(); ++i) {
        Dictionary spawn = existing_spawns[i];
        String type = String(spawn.get("spawn_type", spawn.get("type", spawn.get("name", "")))).to_lower();
        const bool is_player = type == "playerstart" || type == "playerspawn" || type == "player" || type.contains("player");
        const bool is_portal = type == "portal" || type == "exit_portal";
        if (replace_authored_route_markers && ((emit_player_start && is_player) || (emit_exit_portals && is_portal))) continue;
        if (String(spawn.get("source", "")) == "semantic") continue;
        preserved_spawns.append(spawn);
    }

    Array semantic_spawns;
    bool player_emitted = false;
    Array rooms = context.get("rooms", Array());
    for (int i = 0; i < rooms.size(); ++i) {
        Dictionary room = rooms[i];
        String role = String(room.get("spatial_role", "")).to_lower();
        String type = String(room.get("type", "")).to_lower();
        PackedStringArray tags = room.get("gameplay_tags", PackedStringArray());
        bool is_entry = role == "entry" || type.contains("entry") || type.contains("start");
        bool is_exit = role == "exit" || type.contains("exit");
        bool is_encounter = type.contains("encounter") || type.contains("boss");
        for (int tag_index = 0; tag_index < tags.size(); ++tag_index) {
            String tag = tags[tag_index].to_lower();
            is_entry = is_entry || tag == "entry" || tag == "player_start" || tag == "start";
            is_exit = is_exit || tag == "exit" || tag == "portal";
            is_encounter = is_encounter || tag == "encounter" || tag == "boss";
        }

        Vector3 position = room.get("position", Vector3());
        Vector3 size = room.get("size", Vector3());
        Vector3 floor_center = position + Vector3(size.x * 0.5f, floor_offset, size.z * 0.5f);
        if (emit_player_start && is_entry && !player_emitted) {
            Dictionary marker;
            marker["position"] = floor_center;
            marker["spawn_type"] = "PlayerStart";
            marker["source"] = "semantic";
            marker["space_id"] = room.get("id", "");
            semantic_spawns.append(marker);
            player_emitted = true;
        }
        if (emit_exit_portals && is_exit) {
            Dictionary marker;
            marker["position"] = floor_center;
            marker["spawn_type"] = "exit_portal";
            marker["source"] = "semantic";
            marker["space_id"] = room.get("id", "");
            semantic_spawns.append(marker);
        }
        if (emit_encounters && is_encounter) {
            Dictionary marker;
            marker["position"] = floor_center;
            marker["spawn_type"] = "enemy";
            marker["source"] = "semantic";
            marker["space_id"] = room.get("id", "");
            semantic_spawns.append(marker);
        }
    }

    semantic_spawns.append_array(preserved_spawns);
    context["vox_spawns"] = semantic_spawns;
    outputs[0] = context;
}

} // namespace godot
