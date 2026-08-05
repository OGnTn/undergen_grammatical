#include "undergen_spatial_model.h"

#include <godot_cpp/core/math.hpp>
#include <godot_cpp/variant/utility_functions.hpp>
#include <set>

namespace godot {

#define SET_VALUE(field, value) do { field = value; emit_changed(); } while (false)

void UnderGenTopologyNodeData::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_id", "value"), &UnderGenTopologyNodeData::set_id);
    ClassDB::bind_method(D_METHOD("get_id"), &UnderGenTopologyNodeData::get_id);
    ClassDB::bind_method(D_METHOD("set_role", "value"), &UnderGenTopologyNodeData::set_role);
    ClassDB::bind_method(D_METHOD("get_role"), &UnderGenTopologyNodeData::get_role);
    ClassDB::bind_method(D_METHOD("set_gameplay_tags", "value"), &UnderGenTopologyNodeData::set_gameplay_tags);
    ClassDB::bind_method(D_METHOD("get_gameplay_tags"), &UnderGenTopologyNodeData::get_gameplay_tags);
    ClassDB::bind_method(D_METHOD("set_traversable", "value"), &UnderGenTopologyNodeData::set_traversable);
    ClassDB::bind_method(D_METHOD("get_traversable"), &UnderGenTopologyNodeData::get_traversable);
    ClassDB::bind_method(D_METHOD("set_scale", "value"), &UnderGenTopologyNodeData::set_scale);
    ClassDB::bind_method(D_METHOD("get_scale"), &UnderGenTopologyNodeData::get_scale);
    ClassDB::bind_method(D_METHOD("set_openness", "value"), &UnderGenTopologyNodeData::set_openness);
    ClassDB::bind_method(D_METHOD("get_openness"), &UnderGenTopologyNodeData::get_openness);
    ClassDB::bind_method(D_METHOD("set_verticality", "value"), &UnderGenTopologyNodeData::set_verticality);
    ClassDB::bind_method(D_METHOD("get_verticality"), &UnderGenTopologyNodeData::get_verticality);
    ClassDB::bind_method(D_METHOD("set_enclosure", "value"), &UnderGenTopologyNodeData::set_enclosure);
    ClassDB::bind_method(D_METHOD("get_enclosure"), &UnderGenTopologyNodeData::get_enclosure);
    ClassDB::bind_method(D_METHOD("set_prominence", "value"), &UnderGenTopologyNodeData::set_prominence);
    ClassDB::bind_method(D_METHOD("get_prominence"), &UnderGenTopologyNodeData::get_prominence);
    ClassDB::bind_method(D_METHOD("set_preferred_axis", "value"), &UnderGenTopologyNodeData::set_preferred_axis);
    ClassDB::bind_method(D_METHOD("get_preferred_axis"), &UnderGenTopologyNodeData::get_preferred_axis);
    ClassDB::bind_method(D_METHOD("set_elevation_range", "value"), &UnderGenTopologyNodeData::set_elevation_range);
    ClassDB::bind_method(D_METHOD("get_elevation_range"), &UnderGenTopologyNodeData::get_elevation_range);
    ClassDB::bind_method(D_METHOD("set_elevation_band", "value"), &UnderGenTopologyNodeData::set_elevation_band);
    ClassDB::bind_method(D_METHOD("get_elevation_band"), &UnderGenTopologyNodeData::get_elevation_band);
    ClassDB::bind_method(D_METHOD("set_parameters", "value"), &UnderGenTopologyNodeData::set_parameters);
    ClassDB::bind_method(D_METHOD("get_parameters"), &UnderGenTopologyNodeData::get_parameters);

    ADD_PROPERTY(PropertyInfo(Variant::STRING, "id"), "set_id", "get_id");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "role", PROPERTY_HINT_ENUM,
            "Generic,Entry,Reveal,Primary Void,Secondary Void,Pocket,Bay,Throat,Shaft,Chasm,Overlook,Crossing,Vista Target,Anchor Mass,Occluder,Divider Mass,Boundary Route,Vertical Connector,Exit,Slot,Bowl,Dome,Undercroft,Gallery,Buttress,Spine,Island,Canopy,Screen,Dead End,Loop Return"), "set_role", "get_role");
    ADD_PROPERTY(PropertyInfo(Variant::PACKED_STRING_ARRAY, "gameplay_tags"), "set_gameplay_tags", "get_gameplay_tags");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "traversable"), "set_traversable", "get_traversable");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "scale", PROPERTY_HINT_RANGE, "0.1,10.0,0.05"), "set_scale", "get_scale");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "openness", PROPERTY_HINT_RANGE, "0,1,0.01"), "set_openness", "get_openness");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "verticality", PROPERTY_HINT_RANGE, "0,1,0.01"), "set_verticality", "get_verticality");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "enclosure", PROPERTY_HINT_RANGE, "0,1,0.01"), "set_enclosure", "get_enclosure");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "prominence", PROPERTY_HINT_RANGE, "0,1,0.01"), "set_prominence", "get_prominence");
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3, "preferred_axis"), "set_preferred_axis", "get_preferred_axis");
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR2, "elevation_range"), "set_elevation_range", "get_elevation_range");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "elevation_band", PROPERTY_HINT_RANGE, "-8,8,1"), "set_elevation_band", "get_elevation_band");
    ADD_PROPERTY(PropertyInfo(Variant::DICTIONARY, "parameters"), "set_parameters", "get_parameters");

    BIND_ENUM_CONSTANT(ROLE_GENERIC); BIND_ENUM_CONSTANT(ROLE_ENTRY); BIND_ENUM_CONSTANT(ROLE_REVEAL);
    BIND_ENUM_CONSTANT(ROLE_PRIMARY_VOID); BIND_ENUM_CONSTANT(ROLE_SECONDARY_VOID); BIND_ENUM_CONSTANT(ROLE_POCKET);
    BIND_ENUM_CONSTANT(ROLE_BAY); BIND_ENUM_CONSTANT(ROLE_THROAT); BIND_ENUM_CONSTANT(ROLE_SHAFT);
    BIND_ENUM_CONSTANT(ROLE_CHASM); BIND_ENUM_CONSTANT(ROLE_OVERLOOK); BIND_ENUM_CONSTANT(ROLE_CROSSING);
    BIND_ENUM_CONSTANT(ROLE_VISTA_TARGET); BIND_ENUM_CONSTANT(ROLE_ANCHOR_MASS); BIND_ENUM_CONSTANT(ROLE_OCCLUDER);
    BIND_ENUM_CONSTANT(ROLE_DIVIDER_MASS); BIND_ENUM_CONSTANT(ROLE_BOUNDARY_ROUTE);
    BIND_ENUM_CONSTANT(ROLE_VERTICAL_CONNECTOR); BIND_ENUM_CONSTANT(ROLE_EXIT);
    BIND_ENUM_CONSTANT(ROLE_SLOT); BIND_ENUM_CONSTANT(ROLE_BOWL); BIND_ENUM_CONSTANT(ROLE_DOME);
    BIND_ENUM_CONSTANT(ROLE_UNDERCROFT); BIND_ENUM_CONSTANT(ROLE_GALLERY); BIND_ENUM_CONSTANT(ROLE_BUTTRESS);
    BIND_ENUM_CONSTANT(ROLE_SPINE); BIND_ENUM_CONSTANT(ROLE_ISLAND); BIND_ENUM_CONSTANT(ROLE_CANOPY);
    BIND_ENUM_CONSTANT(ROLE_SCREEN); BIND_ENUM_CONSTANT(ROLE_DEAD_END); BIND_ENUM_CONSTANT(ROLE_LOOP_RETURN);
}

void UnderGenTopologyNodeData::set_id(const String &v) { SET_VALUE(id, v); } String UnderGenTopologyNodeData::get_id() const { return id; }
void UnderGenTopologyNodeData::set_role(SpaceRole v) { SET_VALUE(role, v); } UnderGenTopologyNodeData::SpaceRole UnderGenTopologyNodeData::get_role() const { return role; }
void UnderGenTopologyNodeData::set_gameplay_tags(const PackedStringArray &v) { SET_VALUE(gameplay_tags, v); } PackedStringArray UnderGenTopologyNodeData::get_gameplay_tags() const { return gameplay_tags; }
void UnderGenTopologyNodeData::set_traversable(bool v) { SET_VALUE(traversable, v); } bool UnderGenTopologyNodeData::get_traversable() const { return traversable; }
void UnderGenTopologyNodeData::set_scale(float v) { SET_VALUE(scale, Math::max(v, 0.01f)); } float UnderGenTopologyNodeData::get_scale() const { return scale; }
void UnderGenTopologyNodeData::set_openness(float v) { SET_VALUE(openness, Math::clamp(v, 0.0f, 1.0f)); } float UnderGenTopologyNodeData::get_openness() const { return openness; }
void UnderGenTopologyNodeData::set_verticality(float v) { SET_VALUE(verticality, Math::clamp(v, 0.0f, 1.0f)); } float UnderGenTopologyNodeData::get_verticality() const { return verticality; }
void UnderGenTopologyNodeData::set_enclosure(float v) { SET_VALUE(enclosure, Math::clamp(v, 0.0f, 1.0f)); } float UnderGenTopologyNodeData::get_enclosure() const { return enclosure; }
void UnderGenTopologyNodeData::set_prominence(float v) { SET_VALUE(prominence, Math::clamp(v, 0.0f, 1.0f)); } float UnderGenTopologyNodeData::get_prominence() const { return prominence; }
void UnderGenTopologyNodeData::set_preferred_axis(const Vector3 &v) { SET_VALUE(preferred_axis, v); } Vector3 UnderGenTopologyNodeData::get_preferred_axis() const { return preferred_axis; }
void UnderGenTopologyNodeData::set_elevation_range(const Vector2 &v) { SET_VALUE(elevation_range, v); } Vector2 UnderGenTopologyNodeData::get_elevation_range() const { return elevation_range; }
void UnderGenTopologyNodeData::set_elevation_band(int v) { SET_VALUE(elevation_band, v); } int UnderGenTopologyNodeData::get_elevation_band() const { return elevation_band; }
void UnderGenTopologyNodeData::set_parameters(const Dictionary &v) { SET_VALUE(parameters, v.duplicate()); } Dictionary UnderGenTopologyNodeData::get_parameters() const { return parameters; }

void UnderGenTopologyEdgeData::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_id", "value"), &UnderGenTopologyEdgeData::set_id); ClassDB::bind_method(D_METHOD("get_id"), &UnderGenTopologyEdgeData::get_id);
    ClassDB::bind_method(D_METHOD("set_from_id", "value"), &UnderGenTopologyEdgeData::set_from_id); ClassDB::bind_method(D_METHOD("get_from_id"), &UnderGenTopologyEdgeData::get_from_id);
    ClassDB::bind_method(D_METHOD("set_to_id", "value"), &UnderGenTopologyEdgeData::set_to_id); ClassDB::bind_method(D_METHOD("get_to_id"), &UnderGenTopologyEdgeData::get_to_id);
    ClassDB::bind_method(D_METHOD("set_relation", "value"), &UnderGenTopologyEdgeData::set_relation); ClassDB::bind_method(D_METHOD("get_relation"), &UnderGenTopologyEdgeData::get_relation);
    ClassDB::bind_method(D_METHOD("set_traversal_type", "value"), &UnderGenTopologyEdgeData::set_traversal_type); ClassDB::bind_method(D_METHOD("get_traversal_type"), &UnderGenTopologyEdgeData::get_traversal_type);
    ClassDB::bind_method(D_METHOD("set_destination_visibility", "value"), &UnderGenTopologyEdgeData::set_destination_visibility); ClassDB::bind_method(D_METHOD("get_destination_visibility"), &UnderGenTopologyEdgeData::get_destination_visibility);
    ClassDB::bind_method(D_METHOD("set_required", "value"), &UnderGenTopologyEdgeData::set_required); ClassDB::bind_method(D_METHOD("get_required"), &UnderGenTopologyEdgeData::get_required);
    ClassDB::bind_method(D_METHOD("set_strength", "value"), &UnderGenTopologyEdgeData::set_strength); ClassDB::bind_method(D_METHOD("get_strength"), &UnderGenTopologyEdgeData::get_strength);
    ClassDB::bind_method(D_METHOD("set_width_start", "value"), &UnderGenTopologyEdgeData::set_width_start); ClassDB::bind_method(D_METHOD("get_width_start"), &UnderGenTopologyEdgeData::get_width_start);
    ClassDB::bind_method(D_METHOD("set_width_end", "value"), &UnderGenTopologyEdgeData::set_width_end); ClassDB::bind_method(D_METHOD("get_width_end"), &UnderGenTopologyEdgeData::get_width_end);
    ClassDB::bind_method(D_METHOD("set_enclosure", "value"), &UnderGenTopologyEdgeData::set_enclosure); ClassDB::bind_method(D_METHOD("get_enclosure"), &UnderGenTopologyEdgeData::get_enclosure);
    ClassDB::bind_method(D_METHOD("set_exposure", "value"), &UnderGenTopologyEdgeData::set_exposure); ClassDB::bind_method(D_METHOD("get_exposure"), &UnderGenTopologyEdgeData::get_exposure);
    ClassDB::bind_method(D_METHOD("set_curvature", "value"), &UnderGenTopologyEdgeData::set_curvature); ClassDB::bind_method(D_METHOD("get_curvature"), &UnderGenTopologyEdgeData::get_curvature);
    ClassDB::bind_method(D_METHOD("set_parameters", "value"), &UnderGenTopologyEdgeData::set_parameters); ClassDB::bind_method(D_METHOD("get_parameters"), &UnderGenTopologyEdgeData::get_parameters);

    ADD_PROPERTY(PropertyInfo(Variant::STRING, "id"), "set_id", "get_id");
    ADD_PROPERTY(PropertyInfo(Variant::STRING, "from_id"), "set_from_id", "get_from_id");
    ADD_PROPERTY(PropertyInfo(Variant::STRING, "to_id"), "set_to_id", "get_to_id");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "relation", PROPERTY_HINT_ENUM, "Connects,Contains,Overlooks,Occludes,Wraps,Crosses,Stacks Above,Visible From,Revealed By"), "set_relation", "get_relation");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "traversal_type", PROPERTY_HINT_ENUM, "None,Boundary,Interior,Crossing,Vertical"), "set_traversal_type", "get_traversal_type");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "destination_visibility", PROPERTY_HINT_ENUM, "Hidden,Partial,Continuous"), "set_destination_visibility", "get_destination_visibility");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "required"), "set_required", "get_required");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "strength", PROPERTY_HINT_RANGE, "0,1,0.01"), "set_strength", "get_strength");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "width_start", PROPERTY_HINT_RANGE, "0.5,32,0.1"), "set_width_start", "get_width_start");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "width_end", PROPERTY_HINT_RANGE, "0.5,32,0.1"), "set_width_end", "get_width_end");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "enclosure", PROPERTY_HINT_RANGE, "0,1,0.01"), "set_enclosure", "get_enclosure");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "exposure", PROPERTY_HINT_RANGE, "0,1,0.01"), "set_exposure", "get_exposure");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "curvature", PROPERTY_HINT_RANGE, "0,1,0.01"), "set_curvature", "get_curvature");
    ADD_PROPERTY(PropertyInfo(Variant::DICTIONARY, "parameters"), "set_parameters", "get_parameters");

    BIND_ENUM_CONSTANT(RELATION_CONNECTS); BIND_ENUM_CONSTANT(RELATION_CONTAINS); BIND_ENUM_CONSTANT(RELATION_OVERLOOKS);
    BIND_ENUM_CONSTANT(RELATION_OCCLUDES); BIND_ENUM_CONSTANT(RELATION_WRAPS); BIND_ENUM_CONSTANT(RELATION_CROSSES);
    BIND_ENUM_CONSTANT(RELATION_STACKS_ABOVE); BIND_ENUM_CONSTANT(RELATION_VISIBLE_FROM); BIND_ENUM_CONSTANT(RELATION_REVEALED_BY);
    BIND_ENUM_CONSTANT(TRAVERSAL_NONE); BIND_ENUM_CONSTANT(TRAVERSAL_BOUNDARY); BIND_ENUM_CONSTANT(TRAVERSAL_INTERIOR);
    BIND_ENUM_CONSTANT(TRAVERSAL_CROSSING); BIND_ENUM_CONSTANT(TRAVERSAL_VERTICAL);
    BIND_ENUM_CONSTANT(VISIBILITY_HIDDEN); BIND_ENUM_CONSTANT(VISIBILITY_PARTIAL); BIND_ENUM_CONSTANT(VISIBILITY_CONTINUOUS);
}

void UnderGenTopologyEdgeData::set_id(const String &v) { SET_VALUE(id, v); } String UnderGenTopologyEdgeData::get_id() const { return id; }
void UnderGenTopologyEdgeData::set_from_id(const String &v) { SET_VALUE(from_id, v); } String UnderGenTopologyEdgeData::get_from_id() const { return from_id; }
void UnderGenTopologyEdgeData::set_to_id(const String &v) { SET_VALUE(to_id, v); } String UnderGenTopologyEdgeData::get_to_id() const { return to_id; }
void UnderGenTopologyEdgeData::set_relation(RelationType v) { SET_VALUE(relation, v); } UnderGenTopologyEdgeData::RelationType UnderGenTopologyEdgeData::get_relation() const { return relation; }
void UnderGenTopologyEdgeData::set_traversal_type(TraversalType v) { SET_VALUE(traversal_type, v); } UnderGenTopologyEdgeData::TraversalType UnderGenTopologyEdgeData::get_traversal_type() const { return traversal_type; }
void UnderGenTopologyEdgeData::set_destination_visibility(DestinationVisibility v) { SET_VALUE(destination_visibility, v); } UnderGenTopologyEdgeData::DestinationVisibility UnderGenTopologyEdgeData::get_destination_visibility() const { return destination_visibility; }
void UnderGenTopologyEdgeData::set_required(bool v) { SET_VALUE(required, v); } bool UnderGenTopologyEdgeData::get_required() const { return required; }
void UnderGenTopologyEdgeData::set_strength(float v) { SET_VALUE(strength, Math::clamp(v, 0.0f, 1.0f)); } float UnderGenTopologyEdgeData::get_strength() const { return strength; }
void UnderGenTopologyEdgeData::set_width_start(float v) { SET_VALUE(width_start, Math::max(v, 0.1f)); } float UnderGenTopologyEdgeData::get_width_start() const { return width_start; }
void UnderGenTopologyEdgeData::set_width_end(float v) { SET_VALUE(width_end, Math::max(v, 0.1f)); } float UnderGenTopologyEdgeData::get_width_end() const { return width_end; }
void UnderGenTopologyEdgeData::set_enclosure(float v) { SET_VALUE(enclosure, Math::clamp(v, 0.0f, 1.0f)); } float UnderGenTopologyEdgeData::get_enclosure() const { return enclosure; }
void UnderGenTopologyEdgeData::set_exposure(float v) { SET_VALUE(exposure, Math::clamp(v, 0.0f, 1.0f)); } float UnderGenTopologyEdgeData::get_exposure() const { return exposure; }
void UnderGenTopologyEdgeData::set_curvature(float v) { SET_VALUE(curvature, Math::clamp(v, 0.0f, 1.0f)); } float UnderGenTopologyEdgeData::get_curvature() const { return curvature; }
void UnderGenTopologyEdgeData::set_parameters(const Dictionary &v) { SET_VALUE(parameters, v.duplicate()); } Dictionary UnderGenTopologyEdgeData::get_parameters() const { return parameters; }

void UnderGenSemanticGraph::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_nodes", "value"), &UnderGenSemanticGraph::set_nodes); ClassDB::bind_method(D_METHOD("get_nodes"), &UnderGenSemanticGraph::get_nodes);
    ClassDB::bind_method(D_METHOD("set_edges", "value"), &UnderGenSemanticGraph::set_edges); ClassDB::bind_method(D_METHOD("get_edges"), &UnderGenSemanticGraph::get_edges);
    ClassDB::bind_method(D_METHOD("set_seed", "value"), &UnderGenSemanticGraph::set_seed); ClassDB::bind_method(D_METHOD("get_seed"), &UnderGenSemanticGraph::get_seed);
    ClassDB::bind_method(D_METHOD("set_metadata", "value"), &UnderGenSemanticGraph::set_metadata); ClassDB::bind_method(D_METHOD("get_metadata"), &UnderGenSemanticGraph::get_metadata);
    ClassDB::bind_method(D_METHOD("add_node", "node"), &UnderGenSemanticGraph::add_node);
    ClassDB::bind_method(D_METHOD("add_edge", "edge"), &UnderGenSemanticGraph::add_edge);
    ClassDB::bind_method(D_METHOD("find_node", "id"), &UnderGenSemanticGraph::find_node);
    ClassDB::bind_method(D_METHOD("get_edges_for", "id", "traversal_only"), &UnderGenSemanticGraph::get_edges_for, DEFVAL(false));
    ClassDB::bind_method(D_METHOD("validate_graph"), &UnderGenSemanticGraph::validate_graph);
    ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "nodes"), "set_nodes", "get_nodes");
    ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "edges"), "set_edges", "get_edges");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "seed"), "set_seed", "get_seed");
    ADD_PROPERTY(PropertyInfo(Variant::DICTIONARY, "metadata"), "set_metadata", "get_metadata");
}

void UnderGenSemanticGraph::set_nodes(const TypedArray<UnderGenTopologyNodeData> &v) { SET_VALUE(nodes, v); }
TypedArray<UnderGenTopologyNodeData> UnderGenSemanticGraph::get_nodes() const { return nodes; }
void UnderGenSemanticGraph::set_edges(const TypedArray<UnderGenTopologyEdgeData> &v) { SET_VALUE(edges, v); }
TypedArray<UnderGenTopologyEdgeData> UnderGenSemanticGraph::get_edges() const { return edges; }
void UnderGenSemanticGraph::set_seed(int64_t v) { SET_VALUE(seed, v); } int64_t UnderGenSemanticGraph::get_seed() const { return seed; }
void UnderGenSemanticGraph::set_metadata(const Dictionary &v) { SET_VALUE(metadata, v.duplicate()); } Dictionary UnderGenSemanticGraph::get_metadata() const { return metadata; }
void UnderGenSemanticGraph::add_node(const Ref<UnderGenTopologyNodeData> &v) { if (v.is_valid()) { nodes.append(v); emit_changed(); } }
void UnderGenSemanticGraph::add_edge(const Ref<UnderGenTopologyEdgeData> &v) { if (v.is_valid()) { edges.append(v); emit_changed(); } }
Ref<UnderGenTopologyNodeData> UnderGenSemanticGraph::find_node(const String &p_id) const {
    for (int i = 0; i < nodes.size(); ++i) { Ref<UnderGenTopologyNodeData> n = nodes[i]; if (n.is_valid() && n->get_id() == p_id) return n; }
    return Ref<UnderGenTopologyNodeData>();
}
TypedArray<UnderGenTopologyEdgeData> UnderGenSemanticGraph::get_edges_for(const String &p_id, bool p_traversal_only) const {
    TypedArray<UnderGenTopologyEdgeData> result;
    for (int i = 0; i < edges.size(); ++i) {
        Ref<UnderGenTopologyEdgeData> e = edges[i];
        if (e.is_null() || (p_traversal_only && e->get_traversal_type() == UnderGenTopologyEdgeData::TRAVERSAL_NONE)) continue;
        if (e->get_from_id() == p_id || e->get_to_id() == p_id) result.append(e);
    }
    return result;
}
PackedStringArray UnderGenSemanticGraph::validate_graph() const {
    PackedStringArray errors;
    std::set<String> ids;
    for (int i = 0; i < nodes.size(); ++i) {
        Ref<UnderGenTopologyNodeData> n = nodes[i];
        if (n.is_null()) { errors.append("Null topology node at index " + String::num_int64(i)); continue; }
        if (n->get_id().is_empty()) errors.append("Topology node has an empty id");
        else if (!ids.insert(n->get_id()).second) errors.append("Duplicate topology node id: " + n->get_id());
    }
    for (int i = 0; i < edges.size(); ++i) {
        Ref<UnderGenTopologyEdgeData> e = edges[i];
        if (e.is_null()) { errors.append("Null topology edge at index " + String::num_int64(i)); continue; }
        if (!ids.count(e->get_from_id())) errors.append("Edge " + e->get_id() + " has missing from node: " + e->get_from_id());
        if (!ids.count(e->get_to_id())) errors.append("Edge " + e->get_id() + " has missing to node: " + e->get_to_id());
        if (e->get_from_id() == e->get_to_id()) errors.append("Edge " + e->get_id() + " is a self relation");
    }
    return errors;
}

void UnderGenPathState::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_t", "value"), &UnderGenPathState::set_t); ClassDB::bind_method(D_METHOD("get_t"), &UnderGenPathState::get_t);
    ClassDB::bind_method(D_METHOD("set_state", "value"), &UnderGenPathState::set_state); ClassDB::bind_method(D_METHOD("get_state"), &UnderGenPathState::get_state);
    ClassDB::bind_method(D_METHOD("set_width", "value"), &UnderGenPathState::set_width); ClassDB::bind_method(D_METHOD("get_width"), &UnderGenPathState::get_width);
    ClassDB::bind_method(D_METHOD("set_height", "value"), &UnderGenPathState::set_height); ClassDB::bind_method(D_METHOD("get_height"), &UnderGenPathState::get_height);
    ClassDB::bind_method(D_METHOD("set_floor_flatness", "value"), &UnderGenPathState::set_floor_flatness); ClassDB::bind_method(D_METHOD("get_floor_flatness"), &UnderGenPathState::get_floor_flatness);
    ClassDB::bind_method(D_METHOD("set_exposure", "value"), &UnderGenPathState::set_exposure); ClassDB::bind_method(D_METHOD("get_exposure"), &UnderGenPathState::get_exposure);
    ClassDB::bind_method(D_METHOD("set_lateral_bias", "value"), &UnderGenPathState::set_lateral_bias); ClassDB::bind_method(D_METHOD("get_lateral_bias"), &UnderGenPathState::get_lateral_bias);
    ClassDB::bind_method(D_METHOD("set_bank", "value"), &UnderGenPathState::set_bank); ClassDB::bind_method(D_METHOD("get_bank"), &UnderGenPathState::get_bank);
    ClassDB::bind_method(D_METHOD("set_local_noise_scale", "value"), &UnderGenPathState::set_local_noise_scale); ClassDB::bind_method(D_METHOD("get_local_noise_scale"), &UnderGenPathState::get_local_noise_scale);
    ClassDB::bind_method(D_METHOD("set_left_wall", "value"), &UnderGenPathState::set_left_wall); ClassDB::bind_method(D_METHOD("get_left_wall"), &UnderGenPathState::get_left_wall);
    ClassDB::bind_method(D_METHOD("set_right_wall", "value"), &UnderGenPathState::set_right_wall); ClassDB::bind_method(D_METHOD("get_right_wall"), &UnderGenPathState::get_right_wall);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "t", PROPERTY_HINT_RANGE, "0,1,0.01"), "set_t", "get_t");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "state", PROPERTY_HINT_ENUM, "Enclosed,Compressed,Boundary,Ledge,Balcony,Bridge,Crevice,Partial Opening,Widened Overlook"), "set_state", "get_state");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "width"), "set_width", "get_width"); ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "height"), "set_height", "get_height");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "floor_flatness", PROPERTY_HINT_RANGE, "0,1,0.01"), "set_floor_flatness", "get_floor_flatness");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "exposure", PROPERTY_HINT_RANGE, "0,1,0.01"), "set_exposure", "get_exposure");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "lateral_bias", PROPERTY_HINT_RANGE, "-1,1,0.01"), "set_lateral_bias", "get_lateral_bias");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "bank", PROPERTY_HINT_RANGE, "-90,90,0.5"), "set_bank", "get_bank");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "local_noise_scale"), "set_local_noise_scale", "get_local_noise_scale");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "left_wall"), "set_left_wall", "get_left_wall"); ADD_PROPERTY(PropertyInfo(Variant::BOOL, "right_wall"), "set_right_wall", "get_right_wall");
    BIND_ENUM_CONSTANT(STATE_ENCLOSED); BIND_ENUM_CONSTANT(STATE_COMPRESSED); BIND_ENUM_CONSTANT(STATE_BOUNDARY);
    BIND_ENUM_CONSTANT(STATE_LEDGE); BIND_ENUM_CONSTANT(STATE_BALCONY); BIND_ENUM_CONSTANT(STATE_BRIDGE); BIND_ENUM_CONSTANT(STATE_CREVICE);
    BIND_ENUM_CONSTANT(STATE_PARTIAL_OPENING); BIND_ENUM_CONSTANT(STATE_WIDENED_OVERLOOK);
}

void UnderGenPathState::set_t(float v) { SET_VALUE(t, Math::clamp(v, 0.0f, 1.0f)); } float UnderGenPathState::get_t() const { return t; }
void UnderGenPathState::set_state(SpatialState v) { SET_VALUE(state, v); } UnderGenPathState::SpatialState UnderGenPathState::get_state() const { return state; }
void UnderGenPathState::set_width(float v) { SET_VALUE(width, Math::max(v, 0.1f)); } float UnderGenPathState::get_width() const { return width; }
void UnderGenPathState::set_height(float v) { SET_VALUE(height, Math::max(v, 0.1f)); } float UnderGenPathState::get_height() const { return height; }
void UnderGenPathState::set_floor_flatness(float v) { SET_VALUE(floor_flatness, Math::clamp(v, 0.0f, 1.0f)); } float UnderGenPathState::get_floor_flatness() const { return floor_flatness; }
void UnderGenPathState::set_exposure(float v) { SET_VALUE(exposure, Math::clamp(v, 0.0f, 1.0f)); } float UnderGenPathState::get_exposure() const { return exposure; }
void UnderGenPathState::set_lateral_bias(float v) { SET_VALUE(lateral_bias, Math::clamp(v, -1.0f, 1.0f)); } float UnderGenPathState::get_lateral_bias() const { return lateral_bias; }
void UnderGenPathState::set_bank(float v) { SET_VALUE(bank, v); } float UnderGenPathState::get_bank() const { return bank; }
void UnderGenPathState::set_local_noise_scale(float v) { SET_VALUE(local_noise_scale, Math::max(v, 0.0f)); } float UnderGenPathState::get_local_noise_scale() const { return local_noise_scale; }
void UnderGenPathState::set_left_wall(bool v) { SET_VALUE(left_wall, v); } bool UnderGenPathState::get_left_wall() const { return left_wall; }
void UnderGenPathState::set_right_wall(bool v) { SET_VALUE(right_wall, v); } bool UnderGenPathState::get_right_wall() const { return right_wall; }

void UnderGenEmbeddedSpace::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_id", "value"), &UnderGenEmbeddedSpace::set_id); ClassDB::bind_method(D_METHOD("get_id"), &UnderGenEmbeddedSpace::get_id);
    ClassDB::bind_method(D_METHOD("set_role", "value"), &UnderGenEmbeddedSpace::set_role); ClassDB::bind_method(D_METHOD("get_role"), &UnderGenEmbeddedSpace::get_role);
    ClassDB::bind_method(D_METHOD("set_shape", "value"), &UnderGenEmbeddedSpace::set_shape); ClassDB::bind_method(D_METHOD("get_shape"), &UnderGenEmbeddedSpace::get_shape);
    ClassDB::bind_method(D_METHOD("set_transform", "value"), &UnderGenEmbeddedSpace::set_transform); ClassDB::bind_method(D_METHOD("get_transform"), &UnderGenEmbeddedSpace::get_transform);
    ClassDB::bind_method(D_METHOD("set_position", "value"), &UnderGenEmbeddedSpace::set_position); ClassDB::bind_method(D_METHOD("get_position"), &UnderGenEmbeddedSpace::get_position);
    ClassDB::bind_method(D_METHOD("set_size", "value"), &UnderGenEmbeddedSpace::set_size); ClassDB::bind_method(D_METHOD("get_size"), &UnderGenEmbeddedSpace::get_size);
    ClassDB::bind_method(D_METHOD("set_gameplay_tags", "value"), &UnderGenEmbeddedSpace::set_gameplay_tags); ClassDB::bind_method(D_METHOD("get_gameplay_tags"), &UnderGenEmbeddedSpace::get_gameplay_tags);
    ClassDB::bind_method(D_METHOD("set_traversable", "value"), &UnderGenEmbeddedSpace::set_traversable); ClassDB::bind_method(D_METHOD("get_traversable"), &UnderGenEmbeddedSpace::get_traversable);
    ClassDB::bind_method(D_METHOD("set_parameters", "value"), &UnderGenEmbeddedSpace::set_parameters); ClassDB::bind_method(D_METHOD("get_parameters"), &UnderGenEmbeddedSpace::get_parameters);
    ClassDB::bind_method(D_METHOD("get_bounds"), &UnderGenEmbeddedSpace::get_bounds);
    ADD_PROPERTY(PropertyInfo(Variant::STRING, "id"), "set_id", "get_id");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "role", PROPERTY_HINT_ENUM, "Generic,Entry,Reveal,Primary Void,Secondary Void,Pocket,Bay,Throat,Shaft,Chasm,Overlook,Crossing,Vista Target,Anchor Mass,Occluder,Divider Mass,Boundary Route,Vertical Connector,Exit"), "set_role", "get_role");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "shape", PROPERTY_HINT_ENUM, "Box,Ellipsoid,Capsule,Mass"), "set_shape", "get_shape");
    ADD_PROPERTY(PropertyInfo(Variant::TRANSFORM3D, "transform"), "set_transform", "get_transform");
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3, "size"), "set_size", "get_size");
    ADD_PROPERTY(PropertyInfo(Variant::PACKED_STRING_ARRAY, "gameplay_tags"), "set_gameplay_tags", "get_gameplay_tags");
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "traversable"), "set_traversable", "get_traversable");
    ADD_PROPERTY(PropertyInfo(Variant::DICTIONARY, "parameters"), "set_parameters", "get_parameters");
    BIND_ENUM_CONSTANT(SHAPE_BOX); BIND_ENUM_CONSTANT(SHAPE_ELLIPSOID); BIND_ENUM_CONSTANT(SHAPE_CAPSULE); BIND_ENUM_CONSTANT(SHAPE_MASS);
}

void UnderGenEmbeddedSpace::set_id(const String &v) { SET_VALUE(id, v); } String UnderGenEmbeddedSpace::get_id() const { return id; }
void UnderGenEmbeddedSpace::set_role(UnderGenTopologyNodeData::SpaceRole v) { SET_VALUE(role, v); } UnderGenTopologyNodeData::SpaceRole UnderGenEmbeddedSpace::get_role() const { return role; }
void UnderGenEmbeddedSpace::set_shape(ShapeType v) { SET_VALUE(shape, v); } UnderGenEmbeddedSpace::ShapeType UnderGenEmbeddedSpace::get_shape() const { return shape; }
void UnderGenEmbeddedSpace::set_transform(const Transform3D &v) { SET_VALUE(transform, v); } Transform3D UnderGenEmbeddedSpace::get_transform() const { return transform; }
void UnderGenEmbeddedSpace::set_position(const Vector3 &v) { transform.origin = v; emit_changed(); } Vector3 UnderGenEmbeddedSpace::get_position() const { return transform.origin; }
void UnderGenEmbeddedSpace::set_size(const Vector3 &v) { SET_VALUE(size, v.abs()); } Vector3 UnderGenEmbeddedSpace::get_size() const { return size; }
void UnderGenEmbeddedSpace::set_gameplay_tags(const PackedStringArray &v) { SET_VALUE(gameplay_tags, v); } PackedStringArray UnderGenEmbeddedSpace::get_gameplay_tags() const { return gameplay_tags; }
void UnderGenEmbeddedSpace::set_traversable(bool v) { SET_VALUE(traversable, v); } bool UnderGenEmbeddedSpace::get_traversable() const { return traversable; }
void UnderGenEmbeddedSpace::set_parameters(const Dictionary &v) { SET_VALUE(parameters, v.duplicate()); } Dictionary UnderGenEmbeddedSpace::get_parameters() const { return parameters; }
AABB UnderGenEmbeddedSpace::get_bounds() const { return AABB(transform.origin - size * 0.5f, size); }

void UnderGenSpatialField::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_field_type", "value"), &UnderGenSpatialField::set_field_type); ClassDB::bind_method(D_METHOD("get_field_type"), &UnderGenSpatialField::get_field_type);
    ClassDB::bind_method(D_METHOD("set_bounds", "value"), &UnderGenSpatialField::set_bounds); ClassDB::bind_method(D_METHOD("get_bounds"), &UnderGenSpatialField::get_bounds);
    ClassDB::bind_method(D_METHOD("set_resolution", "value"), &UnderGenSpatialField::set_resolution); ClassDB::bind_method(D_METHOD("get_resolution"), &UnderGenSpatialField::get_resolution);
    ClassDB::bind_method(D_METHOD("set_values", "value"), &UnderGenSpatialField::set_values); ClassDB::bind_method(D_METHOD("get_values"), &UnderGenSpatialField::get_values);
    ClassDB::bind_method(D_METHOD("set_default_value", "value"), &UnderGenSpatialField::set_default_value); ClassDB::bind_method(D_METHOD("get_default_value"), &UnderGenSpatialField::get_default_value);
    ClassDB::bind_method(D_METHOD("initialize", "type", "bounds", "resolution", "default_value"), &UnderGenSpatialField::initialize, DEFVAL(0.0f));
    ClassDB::bind_method(D_METHOD("set_value", "cell", "value"), &UnderGenSpatialField::set_value); ClassDB::bind_method(D_METHOD("get_value", "cell"), &UnderGenSpatialField::get_value); ClassDB::bind_method(D_METHOD("sample", "world_position"), &UnderGenSpatialField::sample);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "field_type", PROPERTY_HINT_ENUM, "Openness,Verticality,Exposure,Enclosure,Occlusion,Prominence,Connectivity Pressure,Surface Suitability"), "set_field_type", "get_field_type");
    ADD_PROPERTY(PropertyInfo(Variant::AABB, "bounds"), "set_bounds", "get_bounds"); ADD_PROPERTY(PropertyInfo(Variant::VECTOR3I, "resolution"), "set_resolution", "get_resolution");
    ADD_PROPERTY(PropertyInfo(Variant::PACKED_FLOAT32_ARRAY, "values"), "set_values", "get_values"); ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "default_value", PROPERTY_HINT_RANGE, "0,1,0.01"), "set_default_value", "get_default_value");
    BIND_ENUM_CONSTANT(FIELD_OPENNESS); BIND_ENUM_CONSTANT(FIELD_VERTICALITY); BIND_ENUM_CONSTANT(FIELD_EXPOSURE); BIND_ENUM_CONSTANT(FIELD_ENCLOSURE); BIND_ENUM_CONSTANT(FIELD_OCCLUSION); BIND_ENUM_CONSTANT(FIELD_PROMINENCE); BIND_ENUM_CONSTANT(FIELD_CONNECTIVITY_PRESSURE); BIND_ENUM_CONSTANT(FIELD_SURFACE_SUITABILITY);
}
void UnderGenSpatialField::set_field_type(FieldType v) { SET_VALUE(field_type, v); } UnderGenSpatialField::FieldType UnderGenSpatialField::get_field_type() const { return field_type; }
void UnderGenSpatialField::set_bounds(const AABB &v) { SET_VALUE(bounds, v); } AABB UnderGenSpatialField::get_bounds() const { return bounds; }
void UnderGenSpatialField::set_resolution(const Vector3i &v) { resolution = v.max(Vector3i(1, 1, 1)); values.resize(resolution.x * resolution.y * resolution.z); values.fill(default_value); emit_changed(); } Vector3i UnderGenSpatialField::get_resolution() const { return resolution; }
void UnderGenSpatialField::set_values(const PackedFloat32Array &v) { SET_VALUE(values, v); } PackedFloat32Array UnderGenSpatialField::get_values() const { return values; }
void UnderGenSpatialField::set_default_value(float v) { SET_VALUE(default_value, Math::clamp(v, 0.0f, 1.0f)); } float UnderGenSpatialField::get_default_value() const { return default_value; }
void UnderGenSpatialField::initialize(FieldType p_type, const AABB &p_bounds, const Vector3i &p_resolution, float p_default) { field_type = p_type; bounds = p_bounds; default_value = Math::clamp(p_default, 0.0f, 1.0f); set_resolution(p_resolution); }
bool UnderGenSpatialField::set_value(const Vector3i &p_cell, float p_value) { if (p_cell.x < 0 || p_cell.y < 0 || p_cell.z < 0 || p_cell.x >= resolution.x || p_cell.y >= resolution.y || p_cell.z >= resolution.z) return false; int index = p_cell.x + resolution.x * (p_cell.y + resolution.y * p_cell.z); if (index >= values.size()) return false; values.set(index, Math::clamp(p_value, 0.0f, 1.0f)); return true; }
float UnderGenSpatialField::get_value(const Vector3i &p_cell) const { if (p_cell.x < 0 || p_cell.y < 0 || p_cell.z < 0 || p_cell.x >= resolution.x || p_cell.y >= resolution.y || p_cell.z >= resolution.z) return default_value; int index = p_cell.x + resolution.x * (p_cell.y + resolution.y * p_cell.z); return index < values.size() ? values[index] : default_value; }
float UnderGenSpatialField::sample(const Vector3 &p_world_position) const { if (bounds.size.x <= 0.0f || bounds.size.y <= 0.0f || bounds.size.z <= 0.0f) return default_value; Vector3 n = (p_world_position - bounds.position) / bounds.size; n.x = Math::clamp(n.x, 0.0f, 1.0f); n.y = Math::clamp(n.y, 0.0f, 1.0f); n.z = Math::clamp(n.z, 0.0f, 1.0f); Vector3 scaled = n * Vector3(resolution - Vector3i(1, 1, 1)); return get_value(Vector3i((int)Math::round(scaled.x), (int)Math::round(scaled.y), (int)Math::round(scaled.z))); }

void UnderGenEmbeddedPath::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_id", "value"), &UnderGenEmbeddedPath::set_id); ClassDB::bind_method(D_METHOD("get_id"), &UnderGenEmbeddedPath::get_id);
    ClassDB::bind_method(D_METHOD("set_from_id", "value"), &UnderGenEmbeddedPath::set_from_id); ClassDB::bind_method(D_METHOD("get_from_id"), &UnderGenEmbeddedPath::get_from_id);
    ClassDB::bind_method(D_METHOD("set_to_id", "value"), &UnderGenEmbeddedPath::set_to_id); ClassDB::bind_method(D_METHOD("get_to_id"), &UnderGenEmbeddedPath::get_to_id);
    ClassDB::bind_method(D_METHOD("set_traversal_type", "value"), &UnderGenEmbeddedPath::set_traversal_type); ClassDB::bind_method(D_METHOD("get_traversal_type"), &UnderGenEmbeddedPath::get_traversal_type);
    ClassDB::bind_method(D_METHOD("set_points", "value"), &UnderGenEmbeddedPath::set_points); ClassDB::bind_method(D_METHOD("get_points"), &UnderGenEmbeddedPath::get_points);
    ClassDB::bind_method(D_METHOD("set_states", "value"), &UnderGenEmbeddedPath::set_states); ClassDB::bind_method(D_METHOD("get_states"), &UnderGenEmbeddedPath::get_states);
    ClassDB::bind_method(D_METHOD("set_parameters", "value"), &UnderGenEmbeddedPath::set_parameters); ClassDB::bind_method(D_METHOD("get_parameters"), &UnderGenEmbeddedPath::get_parameters);
    ClassDB::bind_method(D_METHOD("get_bounds"), &UnderGenEmbeddedPath::get_bounds);
    ADD_PROPERTY(PropertyInfo(Variant::STRING, "id"), "set_id", "get_id"); ADD_PROPERTY(PropertyInfo(Variant::STRING, "from_id"), "set_from_id", "get_from_id"); ADD_PROPERTY(PropertyInfo(Variant::STRING, "to_id"), "set_to_id", "get_to_id");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "traversal_type", PROPERTY_HINT_ENUM, "None,Boundary,Interior,Crossing,Vertical"), "set_traversal_type", "get_traversal_type");
    ADD_PROPERTY(PropertyInfo(Variant::PACKED_VECTOR3_ARRAY, "points"), "set_points", "get_points"); ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "states"), "set_states", "get_states");
    ADD_PROPERTY(PropertyInfo(Variant::DICTIONARY, "parameters"), "set_parameters", "get_parameters");
}

void UnderGenEmbeddedPath::set_id(const String &v) { SET_VALUE(id, v); } String UnderGenEmbeddedPath::get_id() const { return id; }
void UnderGenEmbeddedPath::set_from_id(const String &v) { SET_VALUE(from_id, v); } String UnderGenEmbeddedPath::get_from_id() const { return from_id; }
void UnderGenEmbeddedPath::set_to_id(const String &v) { SET_VALUE(to_id, v); } String UnderGenEmbeddedPath::get_to_id() const { return to_id; }
void UnderGenEmbeddedPath::set_traversal_type(UnderGenTopologyEdgeData::TraversalType v) { SET_VALUE(traversal_type, v); } UnderGenTopologyEdgeData::TraversalType UnderGenEmbeddedPath::get_traversal_type() const { return traversal_type; }
void UnderGenEmbeddedPath::set_points(const PackedVector3Array &v) { SET_VALUE(points, v); } PackedVector3Array UnderGenEmbeddedPath::get_points() const { return points; }
void UnderGenEmbeddedPath::set_states(const TypedArray<UnderGenPathState> &v) { SET_VALUE(states, v); } TypedArray<UnderGenPathState> UnderGenEmbeddedPath::get_states() const { return states; }
void UnderGenEmbeddedPath::set_parameters(const Dictionary &v) { SET_VALUE(parameters, v.duplicate()); } Dictionary UnderGenEmbeddedPath::get_parameters() const { return parameters; }
AABB UnderGenEmbeddedPath::get_bounds() const {
    if (points.is_empty()) return AABB();
    float margin = parameters.get("radius", 3.0f);
    Vector3 minp = points[0], maxp = points[0];
    for (int i = 1; i < points.size(); ++i) { minp = minp.min(points[i]); maxp = maxp.max(points[i]); }
    Vector3 pad(margin, margin, margin);
    return AABB(minp - pad, (maxp - minp) + pad * 2.0f);
}

void UnderGenEmbeddedLayout::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_source_graph", "value"), &UnderGenEmbeddedLayout::set_source_graph); ClassDB::bind_method(D_METHOD("get_source_graph"), &UnderGenEmbeddedLayout::get_source_graph);
    ClassDB::bind_method(D_METHOD("set_spaces", "value"), &UnderGenEmbeddedLayout::set_spaces); ClassDB::bind_method(D_METHOD("get_spaces"), &UnderGenEmbeddedLayout::get_spaces);
    ClassDB::bind_method(D_METHOD("set_paths", "value"), &UnderGenEmbeddedLayout::set_paths); ClassDB::bind_method(D_METHOD("get_paths"), &UnderGenEmbeddedLayout::get_paths);
    ClassDB::bind_method(D_METHOD("set_fields", "value"), &UnderGenEmbeddedLayout::set_fields); ClassDB::bind_method(D_METHOD("get_fields"), &UnderGenEmbeddedLayout::get_fields);
    ClassDB::bind_method(D_METHOD("set_constraint_report", "value"), &UnderGenEmbeddedLayout::set_constraint_report); ClassDB::bind_method(D_METHOD("get_constraint_report"), &UnderGenEmbeddedLayout::get_constraint_report);
    ClassDB::bind_method(D_METHOD("get_revision"), &UnderGenEmbeddedLayout::get_revision);
    ClassDB::bind_method(D_METHOD("add_space", "space"), &UnderGenEmbeddedLayout::add_space); ClassDB::bind_method(D_METHOD("add_path", "path"), &UnderGenEmbeddedLayout::add_path);
    ClassDB::bind_method(D_METHOD("add_field", "field"), &UnderGenEmbeddedLayout::add_field); ClassDB::bind_method(D_METHOD("find_field", "type"), &UnderGenEmbeddedLayout::find_field);
    ClassDB::bind_method(D_METHOD("find_space", "id"), &UnderGenEmbeddedLayout::find_space);
    ClassDB::bind_method(D_METHOD("move_space", "id", "position"), &UnderGenEmbeddedLayout::move_space);
    ClassDB::bind_method(D_METHOD("set_space_elevation", "id", "elevation"), &UnderGenEmbeddedLayout::set_space_elevation);
    ClassDB::bind_method(D_METHOD("prefer_space_elevation", "id", "elevation", "strength"), &UnderGenEmbeddedLayout::prefer_space_elevation, DEFVAL(0.8f));
    ClassDB::bind_method(D_METHOD("move_elevation_band", "band", "delta_y", "include_structural_spaces"), &UnderGenEmbeddedLayout::move_elevation_band, DEFVAL(false));
    ClassDB::bind_method(D_METHOD("validate_layout", "tolerance"), &UnderGenEmbeddedLayout::validate_layout, DEFVAL(0.05f));
    ClassDB::bind_method(D_METHOD("get_dirty_regions"), &UnderGenEmbeddedLayout::get_dirty_regions); ClassDB::bind_method(D_METHOD("clear_dirty_regions"), &UnderGenEmbeddedLayout::clear_dirty_regions);
    ClassDB::bind_method(D_METHOD("get_bounds"), &UnderGenEmbeddedLayout::get_bounds);
    ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "source_graph", PROPERTY_HINT_RESOURCE_TYPE, "UnderGenSemanticGraph"), "set_source_graph", "get_source_graph");
    ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "spaces"), "set_spaces", "get_spaces"); ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "paths"), "set_paths", "get_paths"); ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "fields"), "set_fields", "get_fields");
    ADD_PROPERTY(PropertyInfo(Variant::DICTIONARY, "constraint_report"), "set_constraint_report", "get_constraint_report");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "revision", PROPERTY_HINT_NONE, "", PROPERTY_USAGE_EDITOR | PROPERTY_USAGE_READ_ONLY), "", "get_revision");
    ADD_SIGNAL(MethodInfo("layout_changed", PropertyInfo(Variant::STRING, "source_id"), PropertyInfo(Variant::AABB, "dirty_bounds"), PropertyInfo(Variant::INT, "revision")));
}

void UnderGenEmbeddedLayout::set_source_graph(const Ref<UnderGenSemanticGraph> &v) { SET_VALUE(source_graph, v); }
Ref<UnderGenSemanticGraph> UnderGenEmbeddedLayout::get_source_graph() const { return source_graph; }
void UnderGenEmbeddedLayout::set_spaces(const TypedArray<UnderGenEmbeddedSpace> &v) { spaces = v; revision++; emit_changed(); }
TypedArray<UnderGenEmbeddedSpace> UnderGenEmbeddedLayout::get_spaces() const { return spaces; }
void UnderGenEmbeddedLayout::set_paths(const TypedArray<UnderGenEmbeddedPath> &v) { paths = v; revision++; emit_changed(); }
TypedArray<UnderGenEmbeddedPath> UnderGenEmbeddedLayout::get_paths() const { return paths; }
void UnderGenEmbeddedLayout::set_fields(const TypedArray<UnderGenSpatialField> &v) { fields = v; revision++; emit_changed(); }
TypedArray<UnderGenSpatialField> UnderGenEmbeddedLayout::get_fields() const { return fields; }
void UnderGenEmbeddedLayout::set_constraint_report(const Dictionary &v) { SET_VALUE(constraint_report, v.duplicate()); }
Dictionary UnderGenEmbeddedLayout::get_constraint_report() const { return constraint_report; }
int64_t UnderGenEmbeddedLayout::get_revision() const { return revision; }
void UnderGenEmbeddedLayout::add_space(const Ref<UnderGenEmbeddedSpace> &v) { if (v.is_valid()) { spaces.append(v); revision++; emit_changed(); } }
void UnderGenEmbeddedLayout::add_path(const Ref<UnderGenEmbeddedPath> &v) { if (v.is_valid()) { paths.append(v); revision++; emit_changed(); } }
void UnderGenEmbeddedLayout::add_field(const Ref<UnderGenSpatialField> &v) { if (v.is_valid()) { fields.append(v); revision++; emit_changed(); } }
Ref<UnderGenSpatialField> UnderGenEmbeddedLayout::find_field(UnderGenSpatialField::FieldType p_type) const { for (int i = 0; i < fields.size(); ++i) { Ref<UnderGenSpatialField> field = fields[i]; if (field.is_valid() && field->get_field_type() == p_type) return field; } return Ref<UnderGenSpatialField>(); }
Ref<UnderGenEmbeddedSpace> UnderGenEmbeddedLayout::find_space(const String &p_id) const {
    for (int i = 0; i < spaces.size(); ++i) { Ref<UnderGenEmbeddedSpace> s = spaces[i]; if (s.is_valid() && s->get_id() == p_id) return s; }
    return Ref<UnderGenEmbeddedSpace>();
}
void UnderGenEmbeddedLayout::_mark_dirty(const AABB &p_bounds) { if (p_bounds.size != Vector3()) dirty_regions.append(p_bounds); }
void UnderGenEmbeddedLayout::_update_connected_path_endpoints(const String &p_id, const Vector3 &p_old, const Vector3 &p_new) {
    for (int i = 0; i < paths.size(); ++i) {
        Ref<UnderGenEmbeddedPath> path = paths[i]; if (path.is_null()) continue;
        PackedVector3Array points = path->get_points(); if (points.is_empty()) continue;
        bool changed = false; AABB old_bounds = path->get_bounds();
        if (path->get_from_id() == p_id) { points.set(0, p_new); changed = true; }
        if (path->get_to_id() == p_id) { points.set(points.size() - 1, p_new); changed = true; }
        if (changed) { path->set_points(points); _mark_dirty(old_bounds.merge(path->get_bounds())); }
    }
}
bool UnderGenEmbeddedLayout::move_space(const String &p_id, const Vector3 &p_position) {
    Ref<UnderGenEmbeddedSpace> space = find_space(p_id); if (space.is_null()) return false;
    AABB old_bounds = space->get_bounds(); Vector3 old_position = space->get_position();
    space->set_position(p_position); _update_connected_path_endpoints(p_id, old_position, p_position);
    AABB dirty = old_bounds.merge(space->get_bounds()); _mark_dirty(dirty); revision++; emit_changed(); emit_signal("layout_changed", p_id, dirty, revision); return true;
}
bool UnderGenEmbeddedLayout::set_space_elevation(const String &p_id, float p_elevation) {
    Ref<UnderGenEmbeddedSpace> space = find_space(p_id); if (space.is_null()) return false;
    Vector3 position = space->get_position(); position.y = p_elevation; return move_space(p_id, position);
}
bool UnderGenEmbeddedLayout::prefer_space_elevation(const String &p_id, float p_elevation, float p_strength) {
    Ref<UnderGenEmbeddedSpace> space = find_space(p_id); if (space.is_null()) return false;
    Dictionary params = space->get_parameters(); params["preferred_elevation"] = p_elevation; params["elevation_strength"] = Math::clamp(p_strength, 0.0f, 1.0f); space->set_parameters(params);
    float current = space->get_position().y; return set_space_elevation(p_id, Math::lerp(current, p_elevation, Math::clamp(p_strength, 0.0f, 1.0f)));
}
bool UnderGenEmbeddedLayout::move_elevation_band(int p_band, float p_delta_y, bool p_include_structural_spaces) {
    if (Math::is_zero_approx(p_delta_y)) return false;
    PackedStringArray moved_ids;
    AABB combined_dirty;
    bool has_dirty = false;
    for (int i = 0; i < spaces.size(); ++i) {
        Ref<UnderGenEmbeddedSpace> space = spaces[i];
        if (space.is_null()) continue;
        Dictionary params = space->get_parameters();
        if ((int)params.get("elevation_band", 0) != p_band) continue;
        const auto role = space->get_role();
        const bool structural = role == UnderGenTopologyNodeData::ROLE_PRIMARY_VOID || role == UnderGenTopologyNodeData::ROLE_ANCHOR_MASS || role == UnderGenTopologyNodeData::ROLE_OCCLUDER || role == UnderGenTopologyNodeData::ROLE_DIVIDER_MASS || role == UnderGenTopologyNodeData::ROLE_BUTTRESS || role == UnderGenTopologyNodeData::ROLE_SPINE || role == UnderGenTopologyNodeData::ROLE_ISLAND || role == UnderGenTopologyNodeData::ROLE_CANOPY || role == UnderGenTopologyNodeData::ROLE_SCREEN;
        if (structural && !p_include_structural_spaces) continue;
        AABB old_bounds = space->get_bounds();
        space->set_position(space->get_position() + Vector3(0, p_delta_y, 0));
        AABB dirty = old_bounds.merge(space->get_bounds());
        combined_dirty = has_dirty ? combined_dirty.merge(dirty) : dirty;
        has_dirty = true;
        moved_ids.append(space->get_id());
    }
    if (moved_ids.is_empty()) return false;
    for (int i = 0; i < paths.size(); ++i) {
        Ref<UnderGenEmbeddedPath> path = paths[i];
        if (path.is_null()) continue;
        const bool move_from = moved_ids.has(path->get_from_id());
        const bool move_to = moved_ids.has(path->get_to_id());
        if (!move_from && !move_to) continue;
        PackedVector3Array points = path->get_points();
        if (points.is_empty()) continue;
        AABB old_bounds = path->get_bounds();
        if (move_from && move_to) {
            for (int p = 0; p < points.size(); ++p) points.set(p, points[p] + Vector3(0, p_delta_y, 0));
        } else if (move_from) {
            points.set(0, points[0] + Vector3(0, p_delta_y, 0));
        } else {
            points.set(points.size() - 1, points[points.size() - 1] + Vector3(0, p_delta_y, 0));
        }
        path->set_points(points);
        AABB dirty = old_bounds.merge(path->get_bounds());
        combined_dirty = has_dirty ? combined_dirty.merge(dirty) : dirty;
        has_dirty = true;
    }
    _mark_dirty(combined_dirty);
    revision++;
    emit_changed();
    emit_signal("layout_changed", "band:" + String::num_int64(p_band), combined_dirty, revision);
    return true;
}
Dictionary UnderGenEmbeddedLayout::validate_layout(float p_tolerance) {
    PackedStringArray violations;
    int checked = 0;
    if (source_graph.is_valid()) {
        TypedArray<UnderGenTopologyEdgeData> graph_edges = source_graph->get_edges();
        for (int i = 0; i < graph_edges.size(); ++i) {
            Ref<UnderGenTopologyEdgeData> edge = graph_edges[i];
            if (edge.is_null()) continue;
            Ref<UnderGenEmbeddedSpace> from = find_space(edge->get_from_id());
            Ref<UnderGenEmbeddedSpace> to = find_space(edge->get_to_id());
            checked++;
            if (from.is_null() || to.is_null()) {
                violations.append("Missing embedded endpoint for " + edge->get_id());
                continue;
            }
            if (edge->get_required() && edge->get_traversal_type() != UnderGenTopologyEdgeData::TRAVERSAL_NONE) {
                bool found_path = false;
                for (int p = 0; p < paths.size(); ++p) {
                    Ref<UnderGenEmbeddedPath> path = paths[p];
                    if (path.is_valid() && path->get_id() == edge->get_id() && path->get_points().size() >= 2) { found_path = true; break; }
                }
                if (!found_path) violations.append("Required traversal has no path: " + edge->get_id());
            }
            if (edge->get_relation() == UnderGenTopologyEdgeData::RELATION_CONTAINS) {
                AABB outer = from->get_bounds().grow(p_tolerance);
                AABB inner = to->get_bounds();
                if (!outer.encloses(inner)) violations.append("Containment violated: " + edge->get_id());
            } else if ((edge->get_relation() == UnderGenTopologyEdgeData::RELATION_STACKS_ABOVE || edge->get_relation() == UnderGenTopologyEdgeData::RELATION_OVERLOOKS) && from->get_position().y <= to->get_position().y + p_tolerance) {
                violations.append("Vertical ordering violated: " + edge->get_id());
            }
        }
    }
    Dictionary report = constraint_report.duplicate();
    const int violation_count = (int)violations.size();
    report["checked_constraints"] = checked;
    report["violations"] = violations;
    report["satisfied_constraints"] = Math::max(0, checked - violation_count);
    report["score"] = checked == 0 ? 1.0f : (float)Math::max(0, checked - violation_count) / (float)checked;
    report["revision"] = revision;
    report["valid"] = violations.is_empty();
    constraint_report = report;
    emit_changed();
    return report;
}
Array UnderGenEmbeddedLayout::get_dirty_regions() const { return dirty_regions; }
void UnderGenEmbeddedLayout::clear_dirty_regions() { dirty_regions.clear(); }
AABB UnderGenEmbeddedLayout::get_bounds() const {
    AABB result; bool initialized = false;
    for (int i = 0; i < spaces.size(); ++i) { Ref<UnderGenEmbeddedSpace> s = spaces[i]; if (s.is_null()) continue; if (!initialized) { result = s->get_bounds(); initialized = true; } else result = result.merge(s->get_bounds()); }
    for (int i = 0; i < paths.size(); ++i) { Ref<UnderGenEmbeddedPath> p = paths[i]; if (p.is_null()) continue; if (!initialized) { result = p->get_bounds(); initialized = true; } else result = result.merge(p->get_bounds()); }
    return result;
}

void UnderGenGeometryOperation::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_id", "value"), &UnderGenGeometryOperation::set_id); ClassDB::bind_method(D_METHOD("get_id"), &UnderGenGeometryOperation::get_id);
    ClassDB::bind_method(D_METHOD("set_source_id", "value"), &UnderGenGeometryOperation::set_source_id); ClassDB::bind_method(D_METHOD("get_source_id"), &UnderGenGeometryOperation::get_source_id);
    ClassDB::bind_method(D_METHOD("set_operation_type", "value"), &UnderGenGeometryOperation::set_operation_type); ClassDB::bind_method(D_METHOD("get_operation_type"), &UnderGenGeometryOperation::get_operation_type);
    ClassDB::bind_method(D_METHOD("set_primitive_type", "value"), &UnderGenGeometryOperation::set_primitive_type); ClassDB::bind_method(D_METHOD("get_primitive_type"), &UnderGenGeometryOperation::get_primitive_type);
    ClassDB::bind_method(D_METHOD("set_transform", "value"), &UnderGenGeometryOperation::set_transform); ClassDB::bind_method(D_METHOD("get_transform"), &UnderGenGeometryOperation::get_transform);
    ClassDB::bind_method(D_METHOD("set_size", "value"), &UnderGenGeometryOperation::set_size); ClassDB::bind_method(D_METHOD("get_size"), &UnderGenGeometryOperation::get_size);
    ClassDB::bind_method(D_METHOD("set_points", "value"), &UnderGenGeometryOperation::set_points); ClassDB::bind_method(D_METHOD("get_points"), &UnderGenGeometryOperation::get_points);
    ClassDB::bind_method(D_METHOD("set_radius", "value"), &UnderGenGeometryOperation::set_radius); ClassDB::bind_method(D_METHOD("get_radius"), &UnderGenGeometryOperation::get_radius);
    ClassDB::bind_method(D_METHOD("set_height", "value"), &UnderGenGeometryOperation::set_height); ClassDB::bind_method(D_METHOD("get_height"), &UnderGenGeometryOperation::get_height);
    ClassDB::bind_method(D_METHOD("set_material_id", "value"), &UnderGenGeometryOperation::set_material_id); ClassDB::bind_method(D_METHOD("get_material_id"), &UnderGenGeometryOperation::get_material_id);
    ClassDB::bind_method(D_METHOD("set_zone_name", "value"), &UnderGenGeometryOperation::set_zone_name); ClassDB::bind_method(D_METHOD("get_zone_name"), &UnderGenGeometryOperation::get_zone_name);
    ClassDB::bind_method(D_METHOD("set_parameters", "value"), &UnderGenGeometryOperation::set_parameters); ClassDB::bind_method(D_METHOD("get_parameters"), &UnderGenGeometryOperation::get_parameters);
    ClassDB::bind_method(D_METHOD("get_bounds"), &UnderGenGeometryOperation::get_bounds);
    ADD_PROPERTY(PropertyInfo(Variant::STRING, "id"), "set_id", "get_id"); ADD_PROPERTY(PropertyInfo(Variant::STRING, "source_id"), "set_source_id", "get_source_id");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "operation_type", PROPERTY_HINT_ENUM, "Subtract Void,Add Mass,Route Clearance,Ledge,Bridge,Undercut,Terrace"), "set_operation_type", "get_operation_type");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "primitive_type", PROPERTY_HINT_ENUM, "Box,Ellipsoid,Capsule,Sweep"), "set_primitive_type", "get_primitive_type");
    ADD_PROPERTY(PropertyInfo(Variant::TRANSFORM3D, "transform"), "set_transform", "get_transform"); ADD_PROPERTY(PropertyInfo(Variant::VECTOR3, "size"), "set_size", "get_size");
    ADD_PROPERTY(PropertyInfo(Variant::PACKED_VECTOR3_ARRAY, "points"), "set_points", "get_points"); ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "radius"), "set_radius", "get_radius"); ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "height"), "set_height", "get_height");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "material_id"), "set_material_id", "get_material_id"); ADD_PROPERTY(PropertyInfo(Variant::STRING, "zone_name"), "set_zone_name", "get_zone_name"); ADD_PROPERTY(PropertyInfo(Variant::DICTIONARY, "parameters"), "set_parameters", "get_parameters");
    BIND_ENUM_CONSTANT(OP_SUBTRACT_VOID); BIND_ENUM_CONSTANT(OP_ADD_MASS); BIND_ENUM_CONSTANT(OP_ROUTE_CLEARANCE); BIND_ENUM_CONSTANT(OP_LEDGE); BIND_ENUM_CONSTANT(OP_BRIDGE); BIND_ENUM_CONSTANT(OP_UNDERCUT); BIND_ENUM_CONSTANT(OP_TERRACE);
    BIND_ENUM_CONSTANT(PRIMITIVE_BOX); BIND_ENUM_CONSTANT(PRIMITIVE_ELLIPSOID); BIND_ENUM_CONSTANT(PRIMITIVE_CAPSULE); BIND_ENUM_CONSTANT(PRIMITIVE_SWEEP);
}

void UnderGenGeometryOperation::set_id(const String &v) { SET_VALUE(id, v); } String UnderGenGeometryOperation::get_id() const { return id; }
void UnderGenGeometryOperation::set_source_id(const String &v) { SET_VALUE(source_id, v); } String UnderGenGeometryOperation::get_source_id() const { return source_id; }
void UnderGenGeometryOperation::set_operation_type(OperationType v) { SET_VALUE(operation_type, v); } UnderGenGeometryOperation::OperationType UnderGenGeometryOperation::get_operation_type() const { return operation_type; }
void UnderGenGeometryOperation::set_primitive_type(PrimitiveType v) { SET_VALUE(primitive_type, v); } UnderGenGeometryOperation::PrimitiveType UnderGenGeometryOperation::get_primitive_type() const { return primitive_type; }
void UnderGenGeometryOperation::set_transform(const Transform3D &v) { SET_VALUE(transform, v); } Transform3D UnderGenGeometryOperation::get_transform() const { return transform; }
void UnderGenGeometryOperation::set_size(const Vector3 &v) { SET_VALUE(size, v.abs()); } Vector3 UnderGenGeometryOperation::get_size() const { return size; }
void UnderGenGeometryOperation::set_points(const PackedVector3Array &v) { SET_VALUE(points, v); } PackedVector3Array UnderGenGeometryOperation::get_points() const { return points; }
void UnderGenGeometryOperation::set_radius(float v) { SET_VALUE(radius, Math::max(v, 0.1f)); } float UnderGenGeometryOperation::get_radius() const { return radius; }
void UnderGenGeometryOperation::set_height(float v) { SET_VALUE(height, Math::max(v, 0.1f)); } float UnderGenGeometryOperation::get_height() const { return height; }
void UnderGenGeometryOperation::set_material_id(int v) { SET_VALUE(material_id, v); } int UnderGenGeometryOperation::get_material_id() const { return material_id; }
void UnderGenGeometryOperation::set_zone_name(const String &v) { SET_VALUE(zone_name, v); } String UnderGenGeometryOperation::get_zone_name() const { return zone_name; }
void UnderGenGeometryOperation::set_parameters(const Dictionary &v) { SET_VALUE(parameters, v.duplicate()); } Dictionary UnderGenGeometryOperation::get_parameters() const { return parameters; }
AABB UnderGenGeometryOperation::get_bounds() const {
    if (primitive_type == PRIMITIVE_SWEEP && !points.is_empty()) {
        Vector3 minp = points[0], maxp = points[0]; for (int i = 1; i < points.size(); ++i) { minp = minp.min(points[i]); maxp = maxp.max(points[i]); }
        float lateral_radius = radius;
        if (!(bool)parameters.get("left_wall", true) || !(bool)parameters.get("right_wall", true)) lateral_radius *= 1.35f;
        Vector3 pad(lateral_radius, Math::max(lateral_radius, height), lateral_radius); return AABB(minp - pad, maxp - minp + pad * 2.0f);
    }
    return AABB(transform.origin - size * 0.5f, size);
}

void UnderGenGeometryPlan::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_source_layout", "value"), &UnderGenGeometryPlan::set_source_layout); ClassDB::bind_method(D_METHOD("get_source_layout"), &UnderGenGeometryPlan::get_source_layout);
    ClassDB::bind_method(D_METHOD("set_operations", "value"), &UnderGenGeometryPlan::set_operations); ClassDB::bind_method(D_METHOD("get_operations"), &UnderGenGeometryPlan::get_operations);
    ClassDB::bind_method(D_METHOD("set_source_revision", "value"), &UnderGenGeometryPlan::set_source_revision); ClassDB::bind_method(D_METHOD("get_source_revision"), &UnderGenGeometryPlan::get_source_revision);
    ClassDB::bind_method(D_METHOD("set_metadata", "value"), &UnderGenGeometryPlan::set_metadata); ClassDB::bind_method(D_METHOD("get_metadata"), &UnderGenGeometryPlan::get_metadata);
    ClassDB::bind_method(D_METHOD("add_operation", "operation"), &UnderGenGeometryPlan::add_operation); ClassDB::bind_method(D_METHOD("remove_operations_for_source", "source_id"), &UnderGenGeometryPlan::remove_operations_for_source); ClassDB::bind_method(D_METHOD("get_bounds"), &UnderGenGeometryPlan::get_bounds);
    ADD_PROPERTY(PropertyInfo(Variant::OBJECT, "source_layout", PROPERTY_HINT_RESOURCE_TYPE, "UnderGenEmbeddedLayout"), "set_source_layout", "get_source_layout"); ADD_PROPERTY(PropertyInfo(Variant::ARRAY, "operations"), "set_operations", "get_operations"); ADD_PROPERTY(PropertyInfo(Variant::INT, "source_revision"), "set_source_revision", "get_source_revision"); ADD_PROPERTY(PropertyInfo(Variant::DICTIONARY, "metadata"), "set_metadata", "get_metadata");
}
void UnderGenGeometryPlan::set_source_layout(const Ref<UnderGenEmbeddedLayout> &v) { SET_VALUE(source_layout, v); } Ref<UnderGenEmbeddedLayout> UnderGenGeometryPlan::get_source_layout() const { return source_layout; }
void UnderGenGeometryPlan::set_operations(const TypedArray<UnderGenGeometryOperation> &v) { SET_VALUE(operations, v); } TypedArray<UnderGenGeometryOperation> UnderGenGeometryPlan::get_operations() const { return operations; }
void UnderGenGeometryPlan::set_source_revision(int64_t v) { SET_VALUE(source_revision, v); } int64_t UnderGenGeometryPlan::get_source_revision() const { return source_revision; }
void UnderGenGeometryPlan::set_metadata(const Dictionary &v) { SET_VALUE(metadata, v.duplicate()); } Dictionary UnderGenGeometryPlan::get_metadata() const { return metadata; }
void UnderGenGeometryPlan::add_operation(const Ref<UnderGenGeometryOperation> &v) { if (v.is_valid()) { operations.append(v); emit_changed(); } }
int UnderGenGeometryPlan::remove_operations_for_source(const String &p_source_id) { int removed = 0; for (int i = operations.size() - 1; i >= 0; --i) { Ref<UnderGenGeometryOperation> op = operations[i]; if (op.is_valid() && op->get_source_id() == p_source_id) { operations.remove_at(i); removed++; } } if (removed) emit_changed(); return removed; }
AABB UnderGenGeometryPlan::get_bounds() const { AABB result; bool initialized = false; for (int i = 0; i < operations.size(); ++i) { Ref<UnderGenGeometryOperation> op = operations[i]; if (op.is_null()) continue; if (!initialized) { result = op->get_bounds(); initialized = true; } else result = result.merge(op->get_bounds()); } return result; }

void UnderGenGenerationPreset::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_preset_name", "value"), &UnderGenGenerationPreset::set_preset_name); ClassDB::bind_method(D_METHOD("get_preset_name"), &UnderGenGenerationPreset::get_preset_name);
    ClassDB::bind_method(D_METHOD("set_builtin_preset", "value"), &UnderGenGenerationPreset::set_builtin_preset); ClassDB::bind_method(D_METHOD("get_builtin_preset"), &UnderGenGenerationPreset::get_builtin_preset);
    ClassDB::bind_method(D_METHOD("set_world_size", "value"), &UnderGenGenerationPreset::set_world_size); ClassDB::bind_method(D_METHOD("get_world_size"), &UnderGenGenerationPreset::get_world_size);
    ClassDB::bind_method(D_METHOD("set_primary_void_size", "value"), &UnderGenGenerationPreset::set_primary_void_size); ClassDB::bind_method(D_METHOD("get_primary_void_size"), &UnderGenGenerationPreset::get_primary_void_size);
    ClassDB::bind_method(D_METHOD("set_anchor_scale", "value"), &UnderGenGenerationPreset::set_anchor_scale); ClassDB::bind_method(D_METHOD("get_anchor_scale"), &UnderGenGenerationPreset::get_anchor_scale);
    ClassDB::bind_method(D_METHOD("set_anchor_eccentricity", "value"), &UnderGenGenerationPreset::set_anchor_eccentricity); ClassDB::bind_method(D_METHOD("get_anchor_eccentricity"), &UnderGenGenerationPreset::get_anchor_eccentricity);
    ClassDB::bind_method(D_METHOD("set_elevation_band_count", "value"), &UnderGenGenerationPreset::set_elevation_band_count); ClassDB::bind_method(D_METHOD("get_elevation_band_count"), &UnderGenGenerationPreset::get_elevation_band_count);
    ClassDB::bind_method(D_METHOD("set_elevation_band_spacing", "value"), &UnderGenGenerationPreset::set_elevation_band_spacing); ClassDB::bind_method(D_METHOD("get_elevation_band_spacing"), &UnderGenGenerationPreset::get_elevation_band_spacing);
    ClassDB::bind_method(D_METHOD("set_secondary_bay_count", "value"), &UnderGenGenerationPreset::set_secondary_bay_count); ClassDB::bind_method(D_METHOD("get_secondary_bay_count"), &UnderGenGenerationPreset::get_secondary_bay_count);
    ClassDB::bind_method(D_METHOD("set_crossing_count", "value"), &UnderGenGenerationPreset::set_crossing_count); ClassDB::bind_method(D_METHOD("get_crossing_count"), &UnderGenGenerationPreset::get_crossing_count);
    ClassDB::bind_method(D_METHOD("set_route_width", "value"), &UnderGenGenerationPreset::set_route_width); ClassDB::bind_method(D_METHOD("get_route_width"), &UnderGenGenerationPreset::get_route_width);
    ClassDB::bind_method(D_METHOD("set_route_height", "value"), &UnderGenGenerationPreset::set_route_height); ClassDB::bind_method(D_METHOD("get_route_height"), &UnderGenGenerationPreset::get_route_height);
    ClassDB::bind_method(D_METHOD("set_openness", "value"), &UnderGenGenerationPreset::set_openness); ClassDB::bind_method(D_METHOD("get_openness"), &UnderGenGenerationPreset::get_openness);
    ClassDB::bind_method(D_METHOD("set_verticality", "value"), &UnderGenGenerationPreset::set_verticality); ClassDB::bind_method(D_METHOD("get_verticality"), &UnderGenGenerationPreset::get_verticality);
    ClassDB::bind_method(D_METHOD("set_exposure", "value"), &UnderGenGenerationPreset::set_exposure); ClassDB::bind_method(D_METHOD("get_exposure"), &UnderGenGenerationPreset::get_exposure);
    ClassDB::bind_method(D_METHOD("set_noise_amplitude", "value"), &UnderGenGenerationPreset::set_noise_amplitude); ClassDB::bind_method(D_METHOD("get_noise_amplitude"), &UnderGenGenerationPreset::get_noise_amplitude);
    ClassDB::bind_method(D_METHOD("configure_builtin", "preset"), &UnderGenGenerationPreset::configure_builtin);
    ADD_PROPERTY(PropertyInfo(Variant::STRING, "preset_name"), "set_preset_name", "get_preset_name"); ADD_PROPERTY(PropertyInfo(Variant::INT, "builtin_preset", PROPERTY_HINT_ENUM, "Layered Chasm,Compact Cave,Branching Cavern"), "set_builtin_preset", "get_builtin_preset");
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3I, "world_size"), "set_world_size", "get_world_size"); ADD_PROPERTY(PropertyInfo(Variant::VECTOR3, "primary_void_size"), "set_primary_void_size", "get_primary_void_size");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "anchor_scale", PROPERTY_HINT_RANGE, "0,0.8,0.01"), "set_anchor_scale", "get_anchor_scale"); ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "anchor_eccentricity", PROPERTY_HINT_RANGE, "-0.8,0.8,0.01"), "set_anchor_eccentricity", "get_anchor_eccentricity");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "elevation_band_count", PROPERTY_HINT_RANGE, "1,8,1"), "set_elevation_band_count", "get_elevation_band_count"); ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "elevation_band_spacing", PROPERTY_HINT_RANGE, "2,40,0.5"), "set_elevation_band_spacing", "get_elevation_band_spacing");
    ADD_PROPERTY(PropertyInfo(Variant::INT, "secondary_bay_count", PROPERTY_HINT_RANGE, "0,12,1"), "set_secondary_bay_count", "get_secondary_bay_count"); ADD_PROPERTY(PropertyInfo(Variant::INT, "crossing_count", PROPERTY_HINT_RANGE, "0,6,1"), "set_crossing_count", "get_crossing_count");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "route_width", PROPERTY_HINT_RANGE, "1,20,0.25"), "set_route_width", "get_route_width"); ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "route_height", PROPERTY_HINT_RANGE, "1,20,0.25"), "set_route_height", "get_route_height");
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "openness", PROPERTY_HINT_RANGE, "0,1,0.01"), "set_openness", "get_openness"); ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "verticality", PROPERTY_HINT_RANGE, "0,1,0.01"), "set_verticality", "get_verticality"); ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "exposure", PROPERTY_HINT_RANGE, "0,1,0.01"), "set_exposure", "get_exposure"); ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "noise_amplitude", PROPERTY_HINT_RANGE, "0,8,0.1"), "set_noise_amplitude", "get_noise_amplitude");
    BIND_ENUM_CONSTANT(PRESET_LAYERED_CHASM); BIND_ENUM_CONSTANT(PRESET_COMPACT_CAVE); BIND_ENUM_CONSTANT(PRESET_BRANCHING_CAVERN);
}

void UnderGenGenerationPreset::set_preset_name(const String &v) { SET_VALUE(preset_name, v); } String UnderGenGenerationPreset::get_preset_name() const { return preset_name; }
void UnderGenGenerationPreset::set_builtin_preset(BuiltinPreset v) { SET_VALUE(builtin_preset, v); } UnderGenGenerationPreset::BuiltinPreset UnderGenGenerationPreset::get_builtin_preset() const { return builtin_preset; }
void UnderGenGenerationPreset::set_world_size(const Vector3i &v) { SET_VALUE(world_size, v.max(Vector3i(16, 16, 16))); } Vector3i UnderGenGenerationPreset::get_world_size() const { return world_size; }
void UnderGenGenerationPreset::set_primary_void_size(const Vector3 &v) { SET_VALUE(primary_void_size, v.abs()); } Vector3 UnderGenGenerationPreset::get_primary_void_size() const { return primary_void_size; }
void UnderGenGenerationPreset::set_anchor_scale(float v) { SET_VALUE(anchor_scale, Math::clamp(v, 0.0f, 0.9f)); } float UnderGenGenerationPreset::get_anchor_scale() const { return anchor_scale; }
void UnderGenGenerationPreset::set_anchor_eccentricity(float v) { SET_VALUE(anchor_eccentricity, Math::clamp(v, -0.9f, 0.9f)); } float UnderGenGenerationPreset::get_anchor_eccentricity() const { return anchor_eccentricity; }
void UnderGenGenerationPreset::set_elevation_band_count(int v) { SET_VALUE(elevation_band_count, Math::clamp(v, 1, 8)); } int UnderGenGenerationPreset::get_elevation_band_count() const { return elevation_band_count; }
void UnderGenGenerationPreset::set_elevation_band_spacing(float v) { SET_VALUE(elevation_band_spacing, Math::max(v, 1.0f)); } float UnderGenGenerationPreset::get_elevation_band_spacing() const { return elevation_band_spacing; }
void UnderGenGenerationPreset::set_secondary_bay_count(int v) { SET_VALUE(secondary_bay_count, Math::max(v, 0)); } int UnderGenGenerationPreset::get_secondary_bay_count() const { return secondary_bay_count; }
void UnderGenGenerationPreset::set_crossing_count(int v) { SET_VALUE(crossing_count, Math::max(v, 0)); } int UnderGenGenerationPreset::get_crossing_count() const { return crossing_count; }
void UnderGenGenerationPreset::set_route_width(float v) { SET_VALUE(route_width, Math::max(v, 0.5f)); } float UnderGenGenerationPreset::get_route_width() const { return route_width; }
void UnderGenGenerationPreset::set_route_height(float v) { SET_VALUE(route_height, Math::max(v, 0.5f)); } float UnderGenGenerationPreset::get_route_height() const { return route_height; }
void UnderGenGenerationPreset::set_openness(float v) { SET_VALUE(openness, Math::clamp(v, 0.0f, 1.0f)); } float UnderGenGenerationPreset::get_openness() const { return openness; }
void UnderGenGenerationPreset::set_verticality(float v) { SET_VALUE(verticality, Math::clamp(v, 0.0f, 1.0f)); } float UnderGenGenerationPreset::get_verticality() const { return verticality; }
void UnderGenGenerationPreset::set_exposure(float v) { SET_VALUE(exposure, Math::clamp(v, 0.0f, 1.0f)); } float UnderGenGenerationPreset::get_exposure() const { return exposure; }
void UnderGenGenerationPreset::set_noise_amplitude(float v) { SET_VALUE(noise_amplitude, Math::max(v, 0.0f)); } float UnderGenGenerationPreset::get_noise_amplitude() const { return noise_amplitude; }
void UnderGenGenerationPreset::configure_builtin(BuiltinPreset p) {
    builtin_preset = p;
    if (p == PRESET_COMPACT_CAVE) {
        preset_name = "Compact Cave"; world_size = Vector3i(128, 64, 128); primary_void_size = Vector3(72, 38, 66); anchor_scale = 0.18f; anchor_eccentricity = 0.15f; elevation_band_count = 2; elevation_band_spacing = 10.0f; secondary_bay_count = 2; crossing_count = 1; route_width = 4.0f; route_height = 5.0f; openness = 0.55f; verticality = 0.4f; exposure = 0.45f; noise_amplitude = 1.2f;
    } else if (p == PRESET_BRANCHING_CAVERN) {
        preset_name = "Branching Cavern"; world_size = Vector3i(224, 96, 224); primary_void_size = Vector3(100, 62, 100); anchor_scale = 0.28f; anchor_eccentricity = -0.2f; elevation_band_count = 3; elevation_band_spacing = 14.0f; secondary_bay_count = 6; crossing_count = 2; route_width = 4.5f; route_height = 5.5f; openness = 0.68f; verticality = 0.62f; exposure = 0.6f; noise_amplitude = 2.0f;
    } else {
        preset_name = "Layered Chasm"; world_size = Vector3i(192, 96, 192); primary_void_size = Vector3(118, 68, 92); anchor_scale = 0.34f; anchor_eccentricity = 0.22f; elevation_band_count = 3; elevation_band_spacing = 16.0f; secondary_bay_count = 3; crossing_count = 1; route_width = 5.0f; route_height = 6.0f; openness = 0.8f; verticality = 0.8f; exposure = 0.7f; noise_amplitude = 1.5f;
    }
    emit_changed();
}

#undef SET_VALUE

} // namespace godot
