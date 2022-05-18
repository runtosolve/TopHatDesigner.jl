using TopHatDesigner, CSV, DataFrames


purlin_data = CSV.read("/Users/crismoen/.julia/dev/TopHatDesigner/database/Purlins.csv",
DataFrame);

top_hat_data = CSV.read("/Users/crismoen/.julia/dev/TopHatDesigner/database/TopHats.csv",
DataFrame);

existing_deck_data = CSV.read("/Users/crismoen/.julia/dev/TopHatDesigner/database/Existing_Deck.csv",
DataFrame);

new_deck_data = CSV.read("/Users/crismoen/.julia/dev/TopHatDesigner/database/New_Deck.csv",
DataFrame);

purlin_type_1 = "Z8x2.5 067"

purlin_type_2 = "none"

purlin_spans = (25.0)  #ft

purlin_size_span_assignment = (1)

purlin_laps = ()

purlin_spacing = 4.0  #ft

frame_flange_width = 16.0  #in

roof_slope = 1/12

existing_deck_type = "SSR Ultra-Dek 24 in. 24 ga "

purlin_line, purlin_line_uplift = UI.existing_roof_UI_mapper(purlin_spans, purlin_laps, purlin_spacing, roof_slope, purlin_data, existing_deck_type, existing_deck_data, frame_flange_width, purlin_type_1, purlin_type_2, purlin_size_span_assignment);


# purlin_line.cross_section_data[1].node_geometry

# using Plots
# plot(purlin_line.cross_section_data[1].node_geometry[:, 1], purlin_line.cross_section_data[1].node_geometry[:,2], markershape = :o, aspect_ratio=:equal)

# minimum(purlin_line.cross_section_data[1].node_geometry[:,2])
# maximum(purlin_line.cross_section_data[1].node_geometry[:,2])

# plot(purlin_line.free_flange_cross_section_data[1].node_geometry[:, 1], purlin_line.free_flange_cross_section_data[1].node_geometry[:,2], markershape = :o, aspect_ratio=:equal)


# purlin_line.applied_pressure*1000*144

# purlin_line.failure_location

# purlin_line.failure_limit_state 

# purlin_line_uplift.applied_pressure*1000*144

# using Plots
# plot(purlin_line.model.inputs.z, purlin_line.internal_forces.T)

top_hat_type = "TH3.5 071"
		
new_deck_type = "SSR SuperLok-Dek 16 in. 24 ga "

# design_code = "AISI S100-16 ASD"


# top_hat_section_index = findfirst(==(top_hat_type), top_hat_data.section_name)

# num_purlin_line_segments = length(purlin_line.inputs.segments)

# top_hat_cross_section_dimensions = Vector{NTuple{21, Float64}}(undef, num_purlin_line_segments)

# for i = 1:num_purlin_line_segments  #for now use the same TopHat everywhere, set up so that TopHat can change in the future
    
#     top_hat_cross_section_dimensions[i] = tuple([top_hat_data[top_hat_section_index, :][i] for i=2:22]...)

# end

# # Define the TopHat material properties.
# top_hat_material_properties = [(29500.0, 0.30, 55.0, 70.0)]; #E, ν, Fy, Fu


# deck_index = findfirst(==(existing_deck_type), existing_deck_data.deck_name)

# # Define the TopHat punchout dimensions.  
# top_hat_punch_out_dimensions = (existing_deck_data[deck_index, 7], existing_deck_data[deck_index, 8]); #length, height

#   #Define the new deck details.

# new_deck_index = findfirst(==(new_deck_type), new_deck_data.deck_name)


# if !ismissing(new_deck_data[new_deck_index, 3]) #screw_fastened

#     new_roof_panel_details = ("screw-fastened", new_deck_data[new_deck_index, 2], new_deck_data[new_deck_index, 3], new_deck_data[new_deck_index, 4], new_deck_data[new_deck_index, 5])

# elseif ismissing(new_deck_data[new_deck_index, 3]) #SSR

#     new_roof_panel_details = ("vertical leg standing seam", new_deck_data[new_deck_index, 7], new_deck_data[new_deck_index, 6], 0.0, 0.0)

# end

# new_roof_panel_material_properties = (29500.0, 0.30, 55.0, 70.0)

# Assemble the purlin line model, now with the addition of TopHat Framing.

# top_hat_purlin_line = TopHatDesigner.define(purlin_line.inputs.design_code, purlin_line.inputs.segments, purlin_line.inputs.spacing, purlin_line.inputs.roof_slope, purlin_line.inputs.cross_section_dimensions, top_hat_cross_section_dimensions, top_hat_punch_out_dimensions, purlin_line.inputs.material_properties, top_hat_material_properties, purlin_line.inputs.deck_details, purlin_line.inputs.deck_material_properties, new_roof_panel_details, new_roof_panel_material_properties, purlin_line.inputs.frame_flange_width, purlin_line.inputs.support_locations, purlin_line.inputs.purlin_frame_connections, purlin_line.inputs.bridging_locations)

# top_hat_purlin_line = TopHatDesigner.define(purlin_line.inputs.design_code, purlin_line.inputs.segments, purlin_line.inputs.spacing, purlin_line.inputs.roof_slope, purlin_line.inputs.cross_section_dimensions, top_hat_cross_section_dimensions, top_hat_punch_out_dimensions, purlin_line.inputs.material_properties, top_hat_material_properties, purlin_line.inputs.deck_details, purlin_line.inputs.deck_material_properties, new_roof_panel_details, new_roof_panel_material_properties, purlin_line.inputs.frame_flange_width, purlin_line.inputs.support_locations, purlin_line.inputs.purlin_frame_connections, purlin_line.inputs.bridging_locations)





    # #Create the TopHatDesigner data structure.
    # top_hat_purlin_line = TopHatDesigner.TopHatDesignerObject()

    # #Add TopHatDesigner user inputs to data structure.

    # top_hat_purlin_line.inputs = TopHatDesigner.Inputs(purlin_line.inputs.design_code, purlin_line.inputs.segments, purlin_line.inputs.spacing, purlin_line.inputs.roof_slope, purlin_line.inputs.cross_section_dimensions, top_hat_cross_section_dimensions, top_hat_punch_out_dimensions, purlin_line.inputs.material_properties, top_hat_material_properties, new_roof_panel_details, new_roof_panel_material_properties, purlin_line.inputs.frame_flange_width, purlin_line.inputs.support_locations, purlin_line.inputs.purlin_frame_connections, purlin_line.inputs.bridging_locations)

    # #Define TopHat cross-section data including nodal geometry, cross-section discretization and section properties.
    # n = [2, 2, 6, 2, 6, 2, 2]
    # n_radius = [3, 3, 3, 3, 3, 3]
    # top_hat_purlin_line.top_hat_cross_section_data = TopHatDesigner.define_top_hat_cross_sections(top_hat_purlin_line.inputs.top_hat_cross_section_dimensions, n, n_radius)

    # #Create the PurlinLine data structure.
    # purlin_line = PurlinLine.PurlinLineObject()

    # #Capture PurlinLine inputs.
    # purlin_line.inputs = PurlinLine.Inputs(design_code, segments, spacing, roof_slope, purlin_cross_section_dimensions, purlin_material_properties, deck_details, deck_material_properties, frame_flange_width, support_locations, purlin_frame_connections, bridging_locations)

    # #Define the purlin cross-section discretization and calculate section properties.
    # n = [2, 2, 5, 2, 2]
    # n_radius = [3, 3, 3, 3]
    # top_hat_purlin_line.purlin_cross_section_data = PurlinLine.define_purlin_section(purlin_line.inputs.cross_section_dimensions, n, n_radius)
    # purlin_line.cross_section_data = top_hat_purlin_line.purlin_cross_section_data

    # #Define the purlin free flange cross-section discretization and calculate section properties.
    # n = [2, 2, 2]
    # n_radius = [3, 3]
    # top_hat_purlin_line.free_flange_cross_section_data = PurlinLine.define_purlin_free_flange_section(purlin_line.inputs.cross_section_dimensions, n, n_radius)
    # purlin_line.free_flange_cross_section_data = top_hat_purlin_line.free_flange_cross_section_data

    # #Define the purlin cross-section discretization for calculating plastic properties.
    # n = [4, 4, 5, 4, 4] .* 10
    # n_radius = [4, 4, 4, 4]
    # purlin_plastic_cross_section_data = PurlinLine.define_purlin_section(purlin_line.inputs.cross_section_dimensions, n, n_radius)
   
    # #Define TopHat cross-section discretizaton for calculating plastic properties.
    # n = [4, 4, 6, 4, 6, 4, 4] .* 10
    # n_radius = [4, 4, 4, 4, 4, 4]
    # top_hat_plastic_cross_section_data = TopHatDesigner.define_top_hat_cross_sections(top_hat_purlin_line.inputs.top_hat_cross_section_dimensions, n, n_radius)

    # #Define TopHat + purlin properties.
    # top_hat_purlin_line.top_hat_purlin_cross_section_data = TopHatDesigner.define_top_hat_purlin_cross_sections(top_hat_purlin_line.inputs.purlin_cross_section_dimensions, top_hat_purlin_line.purlin_cross_section_data, top_hat_purlin_line.top_hat_cross_section_data, purlin_plastic_cross_section_data, top_hat_plastic_cross_section_data)
   
    # #Define TopHat + purlin plastic discretization.
    # top_hat_purlin_plastic_cross_section_data = TopHatDesigner.define_top_hat_purlin_cross_sections(top_hat_purlin_line.inputs.purlin_cross_section_dimensions, purlin_plastic_cross_section_data, top_hat_plastic_cross_section_data, purlin_plastic_cross_section_data, top_hat_plastic_cross_section_data)

    # # #Calculate net section properties and add to TopHatDesigner data structure.
    # # top_hat_purlin_line.top_hat_purlin_net_cross_section_data = TopHatDesigner.define_top_hat_purlin_net_section(top_hat_purlin_line.inputs.purlin_cross_section_dimensions, top_hat_purlin_line.top_hat_cross_section_data, top_hat_plastic_cross_section_data, top_hat_purlin_line.purlin_cross_section_data, purlin_plastic_cross_section_data,  top_hat_purlin_line.top_hat_purlin_cross_section_data, top_hat_purlin_plastic_cross_section_data, top_hat_purlin_line.inputs.top_hat_punch_out_dimensions)

    # top_hat_purlin_node_geometry, top_hat_purlin_element_definitions = TopHatDesigner.generate_top_hat_net_section_purlin_geometry(top_hat_purlin_line.top_hat_purlin_cross_section_data[1], top_hat_purlin_line.top_hat_cross_section_data[1], top_hat_purlin_line.purlin_cross_section_data[1], top_hat_purlin_line.inputs.purlin_cross_section_dimensions[1],top_hat_purlin_line.inputs.top_hat_punch_out_dimensions)


    # purlin_cross_section_data = top_hat_purlin_line.purlin_cross_section_data[1]
    # top_hat_cross_section_data = top_hat_purlin_line.top_hat_cross_section_data[1]

    # num_purlin_nodes = size(purlin_cross_section_data.node_geometry)[1]

    # #Find all cross-section nodes that are in the TopHat punchout region.
    # top_hat_node_geometry = top_hat_cross_section_data.node_geometry[:,1:2]
    # hole_index_y = findall(x->x<=top_hat_punch_out_dimensions[2], top_hat_node_geometry[:,2])

    # using Plots
    # plot(top_hat_node_geometry[:,1], top_hat_node_geometry[:,2], markershape = :o)


    # #Find all the nodes to the left and right of the TopHat webs.
    # n_top_hat_web_minus = top_hat_cross_section_data.n[1] + top_hat_cross_section_data.n[2] + top_hat_cross_section_data.n_radius[1] + top_hat_cross_section_data.n_radius[2] + floor(Int, top_hat_cross_section_data.n[3]/2)

    # n_top_hat_web_plus = sum(top_hat_cross_section_data.n[1:4]) + sum(top_hat_cross_section_data.n_radius[1:4]) + floor(Int, top_hat_cross_section_data.n[5]/2)

    # top_hat_web_x_minus_location = top_hat_node_geometry[n_top_hat_web_minus, 1]
    # top_hat_web_x_plus_location = top_hat_node_geometry[n_top_hat_web_plus, 1]
    # hole_index_x = findall(x->(x<=top_hat_web_x_plus_location) & ((x>=top_hat_web_x_minus_location)), top_hat_node_geometry[:,1])

    # #These are the nodes to be removed.
    # hole_index = sort(intersect(hole_index_y, hole_index_x))

    # #Shift the cross-section nodes to match up with the punch out dimensions.
    # h_purlin = purlin_cross_section_dimensions[5]
    # top_hat_purlin_node_geometry = top_hat_purlin_cross_section_data.node_geometry[:,1:2]
    # top_hat_purlin_node_geometry[[hole_index[1] + num_purlin_nodes, hole_index[end] + num_purlin_nodes], 2] .=  h_purlin + top_hat_punch_out_dimensions[2]

    # #Delete the nodes within the punchout region.
    # remove_hole_nodes_index = hole_index[2:(end-1)] .+ num_purlin_nodes
    # top_hat_purlin_node_geometry = top_hat_purlin_node_geometry[setdiff(1:end, remove_hole_nodes_index), :]

    # #Remove elements at punchouts.
    # top_hat_purlin_element_definitions = top_hat_purlin_cross_section_data.element_definitions
    # remove_hole_elements_index = hole_index[1:end-1] .+ num_purlin_nodes .- 1
    # top_hat_purlin_element_definitions = top_hat_purlin_element_definitions[setdiff(1:end, remove_hole_elements_index), :]

    # #Update nodal connectivity.
    # update_index = remove_hole_elements_index[1]
    # num_removed_elements = length(remove_hole_elements_index)

    # top_hat_purlin_element_definitions[update_index:end, 1:2] = top_hat_purlin_element_definitions[update_index:end, 1:2] .- num_removed_elements .+ 1




# top_hat_purlin_line = TopHatDesigner.define(purlin_line.inputs.design_code, purlin_line.inputs.segments, purlin_line.inputs.spacing, purlin_line.inputs.roof_slope, purlin_line.inputs.cross_section_dimensions, top_hat_cross_section_dimensions, top_hat_punch_out_dimensions, purlin_line.inputs.material_properties, top_hat_material_properties, purlin_line.inputs.deck_details, purlin_line.inputs.deck_material_properties, new_roof_panel_details, new_roof_panel_material_properties, purlin_line.inputs.frame_flange_width, purlin_line.inputs.support_locations, purlin_line.inputs.purlin_frame_connections, purlin_line.inputs.bridging_locations)

# Run a test to calculate the expected roof system failure pressure including the TopHat.

#Run a gravity test.
# top_hat_purlin_line_gravity = deepcopy(top_hat_purlin_line)
# top_hat_purlin_line_gravity.loading_direction = "gravity"
# top_hat_purlin_line_gravity = TopHatDesigner.capacity(top_hat_purlin_line_gravity)

# #Run an uplift test.
# top_hat_purlin_line_uplift = deepcopy(top_hat_purlin_line)
# top_hat_purlin_line_uplift.loading_direction = "uplift"
# top_hat_purlin_line_uplift = TopHatDesigner.capacity(top_hat_purlin_line_uplift)

# return top_hat_purlin_line_gravity, top_hat_purlin_line_uplift

# end

# function generate_purlin_geometry(t, xcoords_center, ycoords_center, roof_slope)

# center_nodes = [xcoords_center ycoords_center zeros(Float64, length(xcoords_center))]

# center_nodes_rotated = Geometry.rotate_nodes(center_nodes, rotation_axis = "z", rotation_center = [0.0, 0.0, 0.0], θ=atan(roof_slope))

#     #Calculate surface normals.
# closed_or_open = 1
# unitnormals = SectionProperties.surface_normals(xcoords_center, ycoords_center, closed_or_open)


# nodenormals = SectionProperties.avg_node_normals(unitnormals, closed_or_open)


# xcoords_out, ycoords_out = SectionProperties.xycoords_along_normal(xcoords_center, ycoords_center, nodenormals, t/2)

# out_nodes = [xcoords_out ycoords_out zeros(Float64, length(xcoords_center))]

# out_nodes_rotated = Geometry.rotate_nodes(out_nodes, rotation_axis = "z", rotation_center = [0.0, 0.0, 0.0], θ=atan(roof_slope))

# xcoords_in, ycoords_in = SectionProperties.xycoords_along_normal(xcoords_center, ycoords_center, nodenormals, -t/2)

# in_nodes = [xcoords_in ycoords_in zeros(Float64, length(xcoords_center))]

# in_nodes_rotated = Geometry.rotate_nodes(in_nodes, rotation_axis = "z", rotation_center = [0.0, 0.0, 0.0], θ=atan(roof_slope))







# ####################

top_hat_purlin_line_gravity, top_hat_purlin_line_uplift = TopHatDesigner.UI.retrofit_UI_mapper(purlin_line, top_hat_data, top_hat_type, existing_deck_type, existing_deck_data, new_deck_type, new_deck_data);

top_hat_purlin_line_gravity.applied_pressure*1000*144

top_hat_purlin_line_uplift.applied_pressure*1000*144

purlin_line.failure_limit_state

top_hat_purlin_line_gravity.Β_distortional_gradient_factor

using Plots
plot(top_hat_purlin_line_gravity.model.z, top_hat_purlin_line_gravity.biaxial_bending_demand_to_capacity.action_Mxx)
plot!(top_hat_purlin_line_gravity.model.z, top_hat_purlin_line_gravity.biaxial_bending_demand_to_capacity.action_Myy)
plot!(top_hat_purlin_line_gravity.model.z, top_hat_purlin_line_gravity.biaxial_bending_demand_to_capacity.action_Myy)

