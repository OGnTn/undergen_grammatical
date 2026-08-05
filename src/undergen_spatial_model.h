#ifndef UNDERGEN_SPATIAL_MODEL_H
#define UNDERGEN_SPATIAL_MODEL_H

#include <godot_cpp/classes/resource.hpp>
#include <godot_cpp/variant/aabb.hpp>
#include <godot_cpp/variant/dictionary.hpp>
#include <godot_cpp/variant/packed_string_array.hpp>
#include <godot_cpp/variant/packed_float32_array.hpp>
#include <godot_cpp/variant/packed_vector3_array.hpp>
#include <godot_cpp/variant/transform3d.hpp>
#include <godot_cpp/variant/typed_array.hpp>
#include <godot_cpp/variant/vector2.hpp>
#include <godot_cpp/variant/vector3i.hpp>

namespace godot {

class UnderGenTopologyNodeData : public Resource {
    GDCLASS(UnderGenTopologyNodeData, Resource);

public:
    enum SpaceRole {
        ROLE_GENERIC,
        ROLE_ENTRY,
        ROLE_REVEAL,
        ROLE_PRIMARY_VOID,
        ROLE_SECONDARY_VOID,
        ROLE_POCKET,
        ROLE_BAY,
        ROLE_THROAT,
        ROLE_SHAFT,
        ROLE_CHASM,
        ROLE_OVERLOOK,
        ROLE_CROSSING,
        ROLE_VISTA_TARGET,
        ROLE_ANCHOR_MASS,
        ROLE_OCCLUDER,
        ROLE_DIVIDER_MASS,
        ROLE_BOUNDARY_ROUTE,
        ROLE_VERTICAL_CONNECTOR,
        ROLE_EXIT,
        ROLE_SLOT,
        ROLE_BOWL,
        ROLE_DOME,
        ROLE_UNDERCROFT,
        ROLE_GALLERY,
        ROLE_BUTTRESS,
        ROLE_SPINE,
        ROLE_ISLAND,
        ROLE_CANOPY,
        ROLE_SCREEN,
        ROLE_DEAD_END,
        ROLE_LOOP_RETURN
    };

private:
    String id;
    SpaceRole role = ROLE_GENERIC;
    PackedStringArray gameplay_tags;
    bool traversable = true;
    float scale = 1.0f;
    float openness = 0.5f;
    float verticality = 0.5f;
    float enclosure = 0.5f;
    float prominence = 0.5f;
    Vector3 preferred_axis = Vector3(0, 1, 0);
    Vector2 elevation_range = Vector2(-1000.0f, 1000.0f);
    int elevation_band = 0;
    Dictionary parameters;

protected:
    static void _bind_methods();

public:
    void set_id(const String &p_value);
    String get_id() const;
    void set_role(SpaceRole p_value);
    SpaceRole get_role() const;
    void set_gameplay_tags(const PackedStringArray &p_value);
    PackedStringArray get_gameplay_tags() const;
    void set_traversable(bool p_value);
    bool get_traversable() const;
    void set_scale(float p_value);
    float get_scale() const;
    void set_openness(float p_value);
    float get_openness() const;
    void set_verticality(float p_value);
    float get_verticality() const;
    void set_enclosure(float p_value);
    float get_enclosure() const;
    void set_prominence(float p_value);
    float get_prominence() const;
    void set_preferred_axis(const Vector3 &p_value);
    Vector3 get_preferred_axis() const;
    void set_elevation_range(const Vector2 &p_value);
    Vector2 get_elevation_range() const;
    void set_elevation_band(int p_value);
    int get_elevation_band() const;
    void set_parameters(const Dictionary &p_value);
    Dictionary get_parameters() const;
};

class UnderGenTopologyEdgeData : public Resource {
    GDCLASS(UnderGenTopologyEdgeData, Resource);

public:
    enum RelationType {
        RELATION_CONNECTS,
        RELATION_CONTAINS,
        RELATION_OVERLOOKS,
        RELATION_OCCLUDES,
        RELATION_WRAPS,
        RELATION_CROSSES,
        RELATION_STACKS_ABOVE,
        RELATION_VISIBLE_FROM,
        RELATION_REVEALED_BY
    };

    enum TraversalType {
        TRAVERSAL_NONE,
        TRAVERSAL_BOUNDARY,
        TRAVERSAL_INTERIOR,
        TRAVERSAL_CROSSING,
        TRAVERSAL_VERTICAL
    };

    enum DestinationVisibility {
        VISIBILITY_HIDDEN,
        VISIBILITY_PARTIAL,
        VISIBILITY_CONTINUOUS
    };

private:
    String id;
    String from_id;
    String to_id;
    RelationType relation = RELATION_CONNECTS;
    TraversalType traversal_type = TRAVERSAL_INTERIOR;
    DestinationVisibility destination_visibility = VISIBILITY_PARTIAL;
    bool required = true;
    float strength = 1.0f;
    float width_start = 4.0f;
    float width_end = 4.0f;
    float enclosure = 0.5f;
    float exposure = 0.5f;
    float curvature = 0.5f;
    Dictionary parameters;

protected:
    static void _bind_methods();

public:
    void set_id(const String &p_value);
    String get_id() const;
    void set_from_id(const String &p_value);
    String get_from_id() const;
    void set_to_id(const String &p_value);
    String get_to_id() const;
    void set_relation(RelationType p_value);
    RelationType get_relation() const;
    void set_traversal_type(TraversalType p_value);
    TraversalType get_traversal_type() const;
    void set_destination_visibility(DestinationVisibility p_value);
    DestinationVisibility get_destination_visibility() const;
    void set_required(bool p_value);
    bool get_required() const;
    void set_strength(float p_value);
    float get_strength() const;
    void set_width_start(float p_value);
    float get_width_start() const;
    void set_width_end(float p_value);
    float get_width_end() const;
    void set_enclosure(float p_value);
    float get_enclosure() const;
    void set_exposure(float p_value);
    float get_exposure() const;
    void set_curvature(float p_value);
    float get_curvature() const;
    void set_parameters(const Dictionary &p_value);
    Dictionary get_parameters() const;
};

class UnderGenSemanticGraph : public Resource {
    GDCLASS(UnderGenSemanticGraph, Resource);

private:
    TypedArray<UnderGenTopologyNodeData> nodes;
    TypedArray<UnderGenTopologyEdgeData> edges;
    int64_t seed = 12345;
    Dictionary metadata;

protected:
    static void _bind_methods();

public:
    void set_nodes(const TypedArray<UnderGenTopologyNodeData> &p_value);
    TypedArray<UnderGenTopologyNodeData> get_nodes() const;
    void set_edges(const TypedArray<UnderGenTopologyEdgeData> &p_value);
    TypedArray<UnderGenTopologyEdgeData> get_edges() const;
    void set_seed(int64_t p_value);
    int64_t get_seed() const;
    void set_metadata(const Dictionary &p_value);
    Dictionary get_metadata() const;

    void add_node(const Ref<UnderGenTopologyNodeData> &p_node);
    void add_edge(const Ref<UnderGenTopologyEdgeData> &p_edge);
    Ref<UnderGenTopologyNodeData> find_node(const String &p_id) const;
    TypedArray<UnderGenTopologyEdgeData> get_edges_for(const String &p_id, bool p_traversal_only = false) const;
    PackedStringArray validate_graph() const;
};

class UnderGenPathState : public Resource {
    GDCLASS(UnderGenPathState, Resource);

public:
    enum SpatialState {
        STATE_ENCLOSED,
        STATE_COMPRESSED,
        STATE_BOUNDARY,
        STATE_LEDGE,
        STATE_BALCONY,
        STATE_BRIDGE,
        STATE_CREVICE,
        STATE_PARTIAL_OPENING,
        STATE_WIDENED_OVERLOOK
    };

private:
    float t = 0.0f;
    SpatialState state = STATE_ENCLOSED;
    float width = 4.0f;
    float height = 4.0f;
    float floor_flatness = 0.8f;
    float exposure = 0.0f;
    float lateral_bias = 0.0f;
    float bank = 0.0f;
    float local_noise_scale = 1.0f;
    bool left_wall = true;
    bool right_wall = true;

protected:
    static void _bind_methods();

public:
    void set_t(float p_value);
    float get_t() const;
    void set_state(SpatialState p_value);
    SpatialState get_state() const;
    void set_width(float p_value);
    float get_width() const;
    void set_height(float p_value);
    float get_height() const;
    void set_floor_flatness(float p_value);
    float get_floor_flatness() const;
    void set_exposure(float p_value);
    float get_exposure() const;
    void set_lateral_bias(float p_value);
    float get_lateral_bias() const;
    void set_bank(float p_value);
    float get_bank() const;
    void set_local_noise_scale(float p_value);
    float get_local_noise_scale() const;
    void set_left_wall(bool p_value);
    bool get_left_wall() const;
    void set_right_wall(bool p_value);
    bool get_right_wall() const;
};

class UnderGenEmbeddedSpace : public Resource {
    GDCLASS(UnderGenEmbeddedSpace, Resource);

public:
    enum ShapeType {
        SHAPE_BOX,
        SHAPE_ELLIPSOID,
        SHAPE_CAPSULE,
        SHAPE_MASS
    };

private:
    String id;
    UnderGenTopologyNodeData::SpaceRole role = UnderGenTopologyNodeData::ROLE_GENERIC;
    ShapeType shape = SHAPE_ELLIPSOID;
    Transform3D transform;
    Vector3 size = Vector3(8, 6, 8);
    PackedStringArray gameplay_tags;
    bool traversable = true;
    Dictionary parameters;

protected:
    static void _bind_methods();

public:
    void set_id(const String &p_value);
    String get_id() const;
    void set_role(UnderGenTopologyNodeData::SpaceRole p_value);
    UnderGenTopologyNodeData::SpaceRole get_role() const;
    void set_shape(ShapeType p_value);
    ShapeType get_shape() const;
    void set_transform(const Transform3D &p_value);
    Transform3D get_transform() const;
    void set_position(const Vector3 &p_value);
    Vector3 get_position() const;
    void set_size(const Vector3 &p_value);
    Vector3 get_size() const;
    void set_gameplay_tags(const PackedStringArray &p_value);
    PackedStringArray get_gameplay_tags() const;
    void set_traversable(bool p_value);
    bool get_traversable() const;
    void set_parameters(const Dictionary &p_value);
    Dictionary get_parameters() const;
    AABB get_bounds() const;
};

class UnderGenSpatialField : public Resource {
    GDCLASS(UnderGenSpatialField, Resource);

public:
    enum FieldType {
        FIELD_OPENNESS,
        FIELD_VERTICALITY,
        FIELD_EXPOSURE,
        FIELD_ENCLOSURE,
        FIELD_OCCLUSION,
        FIELD_PROMINENCE,
        FIELD_CONNECTIVITY_PRESSURE,
        FIELD_SURFACE_SUITABILITY
    };

private:
    FieldType field_type = FIELD_OPENNESS;
    AABB bounds;
    Vector3i resolution = Vector3i(16, 8, 16);
    PackedFloat32Array values;
    float default_value = 0.0f;

protected:
    static void _bind_methods();

public:
    void set_field_type(FieldType p_value);
    FieldType get_field_type() const;
    void set_bounds(const AABB &p_value);
    AABB get_bounds() const;
    void set_resolution(const Vector3i &p_value);
    Vector3i get_resolution() const;
    void set_values(const PackedFloat32Array &p_value);
    PackedFloat32Array get_values() const;
    void set_default_value(float p_value);
    float get_default_value() const;
    void initialize(FieldType p_type, const AABB &p_bounds, const Vector3i &p_resolution, float p_default_value = 0.0f);
    bool set_value(const Vector3i &p_cell, float p_value);
    float get_value(const Vector3i &p_cell) const;
    float sample(const Vector3 &p_world_position) const;
};

class UnderGenEmbeddedPath : public Resource {
    GDCLASS(UnderGenEmbeddedPath, Resource);

private:
    String id;
    String from_id;
    String to_id;
    UnderGenTopologyEdgeData::TraversalType traversal_type = UnderGenTopologyEdgeData::TRAVERSAL_INTERIOR;
    PackedVector3Array points;
    TypedArray<UnderGenPathState> states;
    Dictionary parameters;

protected:
    static void _bind_methods();

public:
    void set_id(const String &p_value);
    String get_id() const;
    void set_from_id(const String &p_value);
    String get_from_id() const;
    void set_to_id(const String &p_value);
    String get_to_id() const;
    void set_traversal_type(UnderGenTopologyEdgeData::TraversalType p_value);
    UnderGenTopologyEdgeData::TraversalType get_traversal_type() const;
    void set_points(const PackedVector3Array &p_value);
    PackedVector3Array get_points() const;
    void set_states(const TypedArray<UnderGenPathState> &p_value);
    TypedArray<UnderGenPathState> get_states() const;
    void set_parameters(const Dictionary &p_value);
    Dictionary get_parameters() const;
    AABB get_bounds() const;
};

class UnderGenEmbeddedLayout : public Resource {
    GDCLASS(UnderGenEmbeddedLayout, Resource);

private:
    Ref<UnderGenSemanticGraph> source_graph;
    TypedArray<UnderGenEmbeddedSpace> spaces;
    TypedArray<UnderGenEmbeddedPath> paths;
    TypedArray<UnderGenSpatialField> fields;
    Dictionary constraint_report;
    Array dirty_regions;
    int64_t revision = 0;

    void _mark_dirty(const AABB &p_bounds);
    void _update_connected_path_endpoints(const String &p_space_id, const Vector3 &p_old_position, const Vector3 &p_new_position);

protected:
    static void _bind_methods();

public:
    void set_source_graph(const Ref<UnderGenSemanticGraph> &p_value);
    Ref<UnderGenSemanticGraph> get_source_graph() const;
    void set_spaces(const TypedArray<UnderGenEmbeddedSpace> &p_value);
    TypedArray<UnderGenEmbeddedSpace> get_spaces() const;
    void set_paths(const TypedArray<UnderGenEmbeddedPath> &p_value);
    TypedArray<UnderGenEmbeddedPath> get_paths() const;
    void set_fields(const TypedArray<UnderGenSpatialField> &p_value);
    TypedArray<UnderGenSpatialField> get_fields() const;
    void set_constraint_report(const Dictionary &p_value);
    Dictionary get_constraint_report() const;
    int64_t get_revision() const;

    void add_space(const Ref<UnderGenEmbeddedSpace> &p_space);
    void add_path(const Ref<UnderGenEmbeddedPath> &p_path);
    void add_field(const Ref<UnderGenSpatialField> &p_field);
    Ref<UnderGenSpatialField> find_field(UnderGenSpatialField::FieldType p_type) const;
    Ref<UnderGenEmbeddedSpace> find_space(const String &p_id) const;
    bool move_space(const String &p_id, const Vector3 &p_position);
    bool set_space_elevation(const String &p_id, float p_elevation);
    bool prefer_space_elevation(const String &p_id, float p_elevation, float p_strength = 0.8f);
    bool move_elevation_band(int p_band, float p_delta_y, bool p_include_structural_spaces = false);
    Dictionary validate_layout(float p_tolerance = 0.05f);
    Array get_dirty_regions() const;
    void clear_dirty_regions();
    AABB get_bounds() const;
};

class UnderGenGeometryOperation : public Resource {
    GDCLASS(UnderGenGeometryOperation, Resource);

public:
    enum OperationType {
        OP_SUBTRACT_VOID,
        OP_ADD_MASS,
        OP_ROUTE_CLEARANCE,
        OP_LEDGE,
        OP_BRIDGE,
        OP_UNDERCUT,
        OP_TERRACE
    };

    enum PrimitiveType {
        PRIMITIVE_BOX,
        PRIMITIVE_ELLIPSOID,
        PRIMITIVE_CAPSULE,
        PRIMITIVE_SWEEP
    };

private:
    String id;
    String source_id;
    OperationType operation_type = OP_SUBTRACT_VOID;
    PrimitiveType primitive_type = PRIMITIVE_ELLIPSOID;
    Transform3D transform;
    Vector3 size = Vector3(8, 6, 8);
    PackedVector3Array points;
    float radius = 2.0f;
    float height = 4.0f;
    int material_id = 0;
    String zone_name;
    Dictionary parameters;

protected:
    static void _bind_methods();

public:
    void set_id(const String &p_value);
    String get_id() const;
    void set_source_id(const String &p_value);
    String get_source_id() const;
    void set_operation_type(OperationType p_value);
    OperationType get_operation_type() const;
    void set_primitive_type(PrimitiveType p_value);
    PrimitiveType get_primitive_type() const;
    void set_transform(const Transform3D &p_value);
    Transform3D get_transform() const;
    void set_size(const Vector3 &p_value);
    Vector3 get_size() const;
    void set_points(const PackedVector3Array &p_value);
    PackedVector3Array get_points() const;
    void set_radius(float p_value);
    float get_radius() const;
    void set_height(float p_value);
    float get_height() const;
    void set_material_id(int p_value);
    int get_material_id() const;
    void set_zone_name(const String &p_value);
    String get_zone_name() const;
    void set_parameters(const Dictionary &p_value);
    Dictionary get_parameters() const;
    AABB get_bounds() const;
};

class UnderGenGeometryPlan : public Resource {
    GDCLASS(UnderGenGeometryPlan, Resource);

private:
    Ref<UnderGenEmbeddedLayout> source_layout;
    TypedArray<UnderGenGeometryOperation> operations;
    int64_t source_revision = 0;
    Dictionary metadata;

protected:
    static void _bind_methods();

public:
    void set_source_layout(const Ref<UnderGenEmbeddedLayout> &p_value);
    Ref<UnderGenEmbeddedLayout> get_source_layout() const;
    void set_operations(const TypedArray<UnderGenGeometryOperation> &p_value);
    TypedArray<UnderGenGeometryOperation> get_operations() const;
    void set_source_revision(int64_t p_value);
    int64_t get_source_revision() const;
    void set_metadata(const Dictionary &p_value);
    Dictionary get_metadata() const;
    void add_operation(const Ref<UnderGenGeometryOperation> &p_operation);
    int remove_operations_for_source(const String &p_source_id);
    AABB get_bounds() const;
};

class UnderGenGenerationPreset : public Resource {
    GDCLASS(UnderGenGenerationPreset, Resource);

public:
    enum BuiltinPreset {
        PRESET_LAYERED_CHASM,
        PRESET_COMPACT_CAVE,
        PRESET_BRANCHING_CAVERN
    };

private:
    String preset_name = "Layered Chasm";
    BuiltinPreset builtin_preset = PRESET_LAYERED_CHASM;
    Vector3i world_size = Vector3i(192, 96, 192);
    Vector3 primary_void_size = Vector3(118, 68, 92);
    float anchor_scale = 0.34f;
    float anchor_eccentricity = 0.22f;
    int elevation_band_count = 3;
    float elevation_band_spacing = 16.0f;
    int secondary_bay_count = 3;
    int crossing_count = 1;
    float route_width = 5.0f;
    float route_height = 6.0f;
    float openness = 0.8f;
    float verticality = 0.8f;
    float exposure = 0.7f;
    float noise_amplitude = 1.5f;

protected:
    static void _bind_methods();

public:
    void set_preset_name(const String &p_value);
    String get_preset_name() const;
    void set_builtin_preset(BuiltinPreset p_value);
    BuiltinPreset get_builtin_preset() const;
    void set_world_size(const Vector3i &p_value);
    Vector3i get_world_size() const;
    void set_primary_void_size(const Vector3 &p_value);
    Vector3 get_primary_void_size() const;
    void set_anchor_scale(float p_value);
    float get_anchor_scale() const;
    void set_anchor_eccentricity(float p_value);
    float get_anchor_eccentricity() const;
    void set_elevation_band_count(int p_value);
    int get_elevation_band_count() const;
    void set_elevation_band_spacing(float p_value);
    float get_elevation_band_spacing() const;
    void set_secondary_bay_count(int p_value);
    int get_secondary_bay_count() const;
    void set_crossing_count(int p_value);
    int get_crossing_count() const;
    void set_route_width(float p_value);
    float get_route_width() const;
    void set_route_height(float p_value);
    float get_route_height() const;
    void set_openness(float p_value);
    float get_openness() const;
    void set_verticality(float p_value);
    float get_verticality() const;
    void set_exposure(float p_value);
    float get_exposure() const;
    void set_noise_amplitude(float p_value);
    float get_noise_amplitude() const;
    void configure_builtin(BuiltinPreset p_preset);
};

} // namespace godot

VARIANT_ENUM_CAST(godot::UnderGenTopologyNodeData::SpaceRole);
VARIANT_ENUM_CAST(godot::UnderGenTopologyEdgeData::RelationType);
VARIANT_ENUM_CAST(godot::UnderGenTopologyEdgeData::TraversalType);
VARIANT_ENUM_CAST(godot::UnderGenTopologyEdgeData::DestinationVisibility);
VARIANT_ENUM_CAST(godot::UnderGenPathState::SpatialState);
VARIANT_ENUM_CAST(godot::UnderGenEmbeddedSpace::ShapeType);
VARIANT_ENUM_CAST(godot::UnderGenSpatialField::FieldType);
VARIANT_ENUM_CAST(godot::UnderGenGeometryOperation::OperationType);
VARIANT_ENUM_CAST(godot::UnderGenGeometryOperation::PrimitiveType);
VARIANT_ENUM_CAST(godot::UnderGenGenerationPreset::BuiltinPreset);

#endif // UNDERGEN_SPATIAL_MODEL_H
