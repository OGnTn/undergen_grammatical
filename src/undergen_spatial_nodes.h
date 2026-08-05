#ifndef UNDERGEN_SPATIAL_NODES_H
#define UNDERGEN_SPATIAL_NODES_H

#include "undergen_node.h"
#include "undergen_spatial_model.h"

namespace godot {

class DensityGrid;

class UnderGenTopologyBuilderNode : public UnderGenNode {
    GDCLASS(UnderGenTopologyBuilderNode, UnderGenNode);

private:
    Ref<UnderGenGenerationPreset> preset;
    bool prefer_legacy_input = true;

    Ref<UnderGenSemanticGraph> _build_preset_graph(int64_t p_seed) const;
    Ref<UnderGenSemanticGraph> _convert_legacy_graph(const Dictionary &p_graph, int64_t p_seed) const;

protected:
    static void _bind_methods();

public:
    void set_preset(const Ref<UnderGenGenerationPreset> &p_value);
    Ref<UnderGenGenerationPreset> get_preset() const;
    void set_prefer_legacy_input(bool p_value);
    bool get_prefer_legacy_input() const;
    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
    virtual Dictionary get_pipeline_input_defaults(const Dictionary &global_inputs) const override;
};

class UnderGenSpatialEmbedderNode : public UnderGenNode {
    GDCLASS(UnderGenSpatialEmbedderNode, UnderGenNode);

private:
    int solver_iterations = 24;
    float constraint_tolerance = 0.05f;
    bool preserve_vertical_bands = true;

    Ref<UnderGenEmbeddedLayout> _embed(const Ref<UnderGenSemanticGraph> &p_graph, int64_t p_seed) const;

protected:
    static void _bind_methods();

public:
    void set_solver_iterations(int p_value);
    int get_solver_iterations() const;
    void set_constraint_tolerance(float p_value);
    float get_constraint_tolerance() const;
    void set_preserve_vertical_bands(bool p_value);
    bool get_preserve_vertical_bands() const;
    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
    virtual Dictionary get_pipeline_input_defaults(const Dictionary &global_inputs) const override;
};

class UnderGenSpatialValidatorNode : public UnderGenNode {
    GDCLASS(UnderGenSpatialValidatorNode, UnderGenNode);

private:
    float minimum_positive_mass_ratio = 0.08f;
    float maximum_positive_mass_ratio = 0.50f;
    float maximum_crossing_ratio = 0.30f;
    int minimum_elevation_bands = 2;
    bool reject_invalid_layout = false;

protected:
    static void _bind_methods();

public:
    void set_minimum_positive_mass_ratio(float p_value);
    float get_minimum_positive_mass_ratio() const;
    void set_maximum_positive_mass_ratio(float p_value);
    float get_maximum_positive_mass_ratio() const;
    void set_maximum_crossing_ratio(float p_value);
    float get_maximum_crossing_ratio() const;
    void set_minimum_elevation_bands(int p_value);
    int get_minimum_elevation_bands() const;
    void set_reject_invalid_layout(bool p_value);
    bool get_reject_invalid_layout() const;
    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
};

class UnderGenGeometryPlannerNode : public UnderGenNode {
    GDCLASS(UnderGenGeometryPlannerNode, UnderGenNode);

private:
    bool generate_boundary_ledges = true;
    bool generate_undercuts = true;
    float ledge_thickness = 2.0f;

    Ref<UnderGenGeometryPlan> _build_plan(const Ref<UnderGenEmbeddedLayout> &p_layout) const;

protected:
    static void _bind_methods();

public:
    void set_generate_boundary_ledges(bool p_value);
    bool get_generate_boundary_ledges() const;
    void set_generate_undercuts(bool p_value);
    bool get_generate_undercuts() const;
    void set_ledge_thickness(float p_value);
    float get_ledge_thickness() const;
    Ref<UnderGenGeometryPlan> build_plan(const Ref<UnderGenEmbeddedLayout> &p_layout) const;
    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
};

class UnderGenGeometryRealizerNode : public UnderGenNode {
    GDCLASS(UnderGenGeometryRealizerNode, UnderGenNode);

private:
    float surface_threshold = 0.0f;
    float voxel_size = 1.0f;
    bool retain_spatial_data = true;

    void _apply_operation(const Ref<DensityGrid> &p_grid, const Ref<UnderGenGeometryOperation> &p_operation, bool p_use_clip = false, const AABB &p_clip = AABB()) const;
    void _apply_primitive(const Ref<DensityGrid> &p_grid, const Ref<UnderGenGeometryOperation> &p_operation, bool p_add_mass, bool p_use_clip, const AABB &p_clip) const;
    void _apply_sweep(const Ref<DensityGrid> &p_grid, const Ref<UnderGenGeometryOperation> &p_operation, bool p_add_mass, bool p_use_clip, const AABB &p_clip) const;
    Dictionary _update_context_from_plan(const Dictionary &p_context, const Ref<UnderGenGeometryPlan> &p_plan, const Ref<DensityGrid> &p_grid) const;

protected:
    static void _bind_methods();

public:
    void set_surface_threshold(float p_value);
    float get_surface_threshold() const;
    void set_voxel_size(float p_value);
    float get_voxel_size() const;
    void set_retain_spatial_data(bool p_value);
    bool get_retain_spatial_data() const;
    Dictionary realize_plan(const Ref<UnderGenGeometryPlan> &p_plan) const;
    Dictionary rebuild_dirty_regions(const Dictionary &p_context, const Ref<UnderGenGeometryPlan> &p_plan, const Array &p_dirty_regions) const;
    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
};

class UnderGenGameplayMarkerNode : public UnderGenNode {
    GDCLASS(UnderGenGameplayMarkerNode, UnderGenNode);

private:
    bool emit_player_start = true;
    bool emit_exit_portals = true;
    bool emit_encounters = false;
    bool replace_authored_route_markers = true;
    float floor_offset = 1.0f;

protected:
    static void _bind_methods();

public:
    void set_emit_player_start(bool p_value);
    bool get_emit_player_start() const;
    void set_emit_exit_portals(bool p_value);
    bool get_emit_exit_portals() const;
    void set_emit_encounters(bool p_value);
    bool get_emit_encounters() const;
    void set_replace_authored_route_markers(bool p_value);
    bool get_replace_authored_route_markers() const;
    void set_floor_offset(float p_value);
    float get_floor_offset() const;
    virtual void _execute(const Dictionary &inputs, Dictionary &outputs) override;
};

} // namespace godot

#endif // UNDERGEN_SPATIAL_NODES_H
