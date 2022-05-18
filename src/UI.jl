module UI

using PurlinLine, Plots, LinesCurvesNodes, CrossSection

using ..TopHatDesigner

function define_lap_section_types(purlin_size_span_assignment)

	#There are 3 possible combinations of purlin types at a lap: 1-1, 2-2, or 1-2.

	num_laps = length(purlin_size_span_assignment) - 1

	lap_section_types = Array{String}(undef, num_laps)

	for i=1:num_laps

		purlin_span_1 = purlin_size_span_assignment[i]
		purlin_span_2 = purlin_size_span_assignment[i+1]

		if (purlin_span_1 == purlin_span_2) & (purlin_span_1 == 1)

			lap_section_types[i] = "1-1"

		elseif (purlin_span_1 == purlin_span_2) & (purlin_span_1 == 2)

			lap_section_types[i] = "2-2"

		elseif purlin_span_1 != purlin_span_2

			lap_section_types[i] = "1-2"

		end

	end

	return lap_section_types

end

function define_lap_section_index(lap_section_types)

	num_lap_sections = length(lap_section_types)
	
	num_unique_lap_sections = length(unique(lap_section_types))

	unique_lap_sections = unique(lap_section_types)

	lap_section_index = Array{Int64}(undef, num_lap_sections)

	for i = 1:num_unique_lap_sections

		index = findall(x->x==unique_lap_sections[i], lap_section_types)

		lap_section_index[index] .= 2 + i

	end

	return lap_section_index

end

function define_lap_segments(purlin_laps, purlin_size_span_assignment)

	num_interior_supports = trunc(Int, length(purlin_laps)/2)

	lap_segments = Array{Tuple{Float64, Float64, Int64, Int64}, 1}(undef, num_interior_supports * 2)

	lap_section_types = define_lap_section_types(purlin_size_span_assignment)

	lap_section_index = define_lap_section_index(lap_section_types)

	for i = 1:num_interior_supports

		lap_segments[2*i - 1] = (purlin_laps[2*i - 1]*12.0, (purlin_laps[2*i - 1]*12.0)/6, lap_section_index[i], 1)
		lap_segments[2*i] = (purlin_laps[2*i]*12.0, (purlin_laps[2*i]*12.0)/6, lap_section_index[i], 1)

	end

	return lap_segments

end

function define_purlin_line_segments(span_segments, lap_segments)

	num_spans = length(span_segments)
	
	num_purlin_line_segments = length(span_segments) + length(lap_segments)

	purlin_line_segments = Array{Tuple{Float64, Float64, Int64, Int64}, 1}(undef, num_purlin_line_segments)

	segment_index = 1

	lap_segment_index = 1	

	for i = 1:num_spans

		if i == 1 #first span
	
			purlin_line_segments[i] = span_segments[1]

			if num_spans > 1
				purlin_line_segments[i+1] = lap_segments[1]
	
				lap_segment_index = lap_segment_index + 1
				segment_index = segment_index + 2
			end
	
		elseif (i > 1) & (i != num_spans) #interior span
	
			purlin_line_segments[segment_index] = lap_segments[lap_segment_index]
			purlin_line_segments[segment_index+1] = span_segments[i]
			purlin_line_segments[segment_index+2] = lap_segments[lap_segment_index+1]
	
			lap_segment_index = lap_segment_index + 2
			segment_index = segment_index + 3
	
		elseif i == num_spans  #end span
	
			purlin_line_segments[segment_index] = lap_segments[lap_segment_index]
			purlin_line_segments[segment_index+1] = span_segments[i]
	
		end

	end

	return purlin_line_segments

end
	
function define_purlin_line_cross_section_dimensions(purlin_line_segments, lap_section_types, purlin_data, purlin_type_1, purlin_type_2)

	purlin_line_cross_section_indices = [purlin_line_segments[i][3] for i=1:length(purlin_line_segments)]

	unique_purlin_line_cross_section_indices = sort(unique(purlin_line_cross_section_indices))

	if !isempty(lap_section_types) #for multiple spans only
	
		if isempty(findall(x->x==2, unique_purlin_line_cross_section_indices)) #if there is no second purlin type defined, then add it, it won't be used though
	
			unique_purlin_line_cross_section_indices = [unique_purlin_line_cross_section_indices[1]; 2; unique_purlin_line_cross_section_indices[2:end]]
	
		end

	end

	num_purlin_line_cross_sections = length(unique_purlin_line_cross_section_indices)

	purlin_line_cross_section_dimensions = Vector{Tuple{String, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64}}(undef, num_purlin_line_cross_sections)

	purlin_index_1 = findfirst(==(purlin_type_1), purlin_data.section_name)
	purlin_index_2 = findfirst(==(purlin_type_2), purlin_data.section_name)
	
	purlin_type_indices = [purlin_index_1; purlin_index_2]

	for i = 1:num_purlin_line_cross_sections

		cross_section_index = unique_purlin_line_cross_section_indices[i]

		if cross_section_index <= 2

			if purlin_type_2 == "none" #not used, set cross-section equal to purlin_type_1
	
				cross_section_index = unique_purlin_line_cross_section_indices[1]
	
			end
			
			purlin_line_cross_section_dimensions[i] = tuple([purlin_data[purlin_type_indices[cross_section_index], :][i] for i=2:17]...)

		elseif cross_section_index > 2 #these are the laps

			lap_type = lap_section_types[cross_section_index - 2]

			if lap_type == "1-1"

				cross_section_index = unique_purlin_line_cross_section_indices[1]

				baseline_cross_section_dimensions = deepcopy(purlin_data[purlin_type_indices[cross_section_index], :])

				#multiply base metal thickness by 2
				baseline_cross_section_dimensions[3] = baseline_cross_section_dimensions[3] * 2.0
				
				purlin_line_cross_section_dimensions[i] = tuple([baseline_cross_section_dimensions[i] for i=2:17]...)

			elseif lap_type == "2-2"

				cross_section_index = unique_purlin_line_cross_section_indices[2]

				baseline_cross_section_dimensions = deepcopy(purlin_data[purlin_type_indices[cross_section_index], :])

				#multiply base metal thickness by 2
				baseline_cross_section_dimensions[3] = baseline_cross_section_dimensions[3] * 2.0

				purlin_line_cross_section_dimensions[i] = tuple([baseline_cross_section_dimensions[i] for i=2:17]...)

			elseif lap_type == "1-2"

				cross_section_index = unique_purlin_line_cross_section_indices[1]

				baseline_cross_section_dimensions = deepcopy(purlin_data[purlin_type_indices[cross_section_index], :])

				#add the purlin 1 and purlin 2 thicknesses together
				baseline_cross_section_dimensions[3] = baseline_cross_section_dimensions[3] + purlin_line_cross_section_dimensions[2][2]

				
				purlin_line_cross_section_dimensions[i] = tuple([baseline_cross_section_dimensions[i] for i=2:17]...)

			end

		end

	end

	return purlin_line_cross_section_dimensions

end

function define_span_segments(purlin_spans, purlin_laps, purlin_size_span_assignment)

	num_spans = length(purlin_spans)

	span_segments = Array{Tuple{Float64, Float64, Int64, Int64}, 1}(undef, num_spans)

	lap_index = 1

	for i = 1:num_spans

		if i == 1 #first spans

			if num_spans == 1 #for single span

				segment_length = purlin_spans[i]*12

			else #multiple spans
				
				segment_length = purlin_spans[i]*12 - purlin_laps[lap_index]*12

			end
			
			span_segments[i] = (segment_length, segment_length/18, purlin_size_span_assignment[i], 1)
	
			lap_index = lap_index + 1
	
		elseif (i > 1) & (i != num_spans) #interior span
	
			segment_length = purlin_spans[i]*12 - purlin_laps[lap_index]*12 - purlin_laps[lap_index+1]*12
			span_segments[i] = (segment_length, segment_length/18, purlin_size_span_assignment[i], 1)
	
			lap_index = lap_index + 2
	
		elseif i==num_spans #end span
	
			segment_length = purlin_spans[i]*12 - purlin_laps[lap_index]*12
			span_segments[i] = (segment_length, segment_length/18, purlin_size_span_assignment[i], 1)
	
		end

	end

	return span_segments

end

function existing_roof_UI_mapper(purlin_spans, purlin_laps, purlin_spacing, roof_slope, purlin_data, existing_deck_type, existing_deck_data, frame_flange_width, purlin_type_1, purlin_type_2, purlin_size_span_assignment)

	design_code = "AISI S100-16 ASD"

	span_segments = define_span_segments(purlin_spans, purlin_laps, purlin_size_span_assignment)

	lap_section_types = define_lap_section_types(purlin_size_span_assignment)

	lap_segments = define_lap_segments(purlin_laps, purlin_size_span_assignment)

	purlin_segments = define_purlin_line_segments(span_segments, lap_segments)

	purlin_spacing = purlin_spacing * 12.0

	roof_slope = rad2deg(atan(roof_slope))

	purlin_cross_section_dimensions = define_purlin_line_cross_section_dimensions(purlin_segments, lap_section_types, purlin_data, purlin_type_1, purlin_type_2)

	purlin_material_properties = [(29500.0, 0.30, 50.0, 70.0)];  #E, ν, Fy, Fu

	deck_index = findfirst(==(existing_deck_type), existing_deck_data.deck_name)
	existing_roof_panel_details = ("screw-fastened", existing_deck_data[deck_index, 2], existing_deck_data[deck_index, 3], existing_deck_data[deck_index, 4], existing_deck_data[deck_index, 5])

	existing_roof_panel_material_properties = (29500.0, 0.30, 55.0, 70.0);  #E, ν, Fy, Fu

	support_locations = [0.0; collect(cumsum(purlin_spans .* 12.0))]

	# if purlin_frame_connection == "Clip-mounted"
		
		purlin_frame_connections = "anti-roll clip"
		
	# elseif purlin_frame_connection == "Direct"
		
		# purlin_frame_connections = "bottom flange connection"
		
	# end

	intermediate_bridging_locations = [ ]

	purlin_line = PurlinLine.build(design_code, purlin_segments, purlin_spacing, roof_slope, purlin_cross_section_dimensions, purlin_material_properties, existing_roof_panel_details, existing_roof_panel_material_properties, frame_flange_width, support_locations, purlin_frame_connections, intermediate_bridging_locations)

	#Run a gravity test.
	purlin_line_gravity = deepcopy(purlin_line)
	purlin_line_gravity.loading_direction = "gravity"
	purlin_line_gravity = PurlinLine.test(purlin_line_gravity)

    #Run an uplift test.
	purlin_line_uplift = deepcopy(purlin_line)
	purlin_line_uplift.loading_direction = "uplift"
	purlin_line_uplift = PurlinLine.test(purlin_line_uplift)

	return purlin_line_gravity, purlin_line_uplift

end

function retrofit_UI_mapper(purlin_line, top_hat_data, top_hat_type, existing_deck_type, existing_deck_data, new_deck_type, new_deck_data)


	top_hat_section_index = findfirst(==(top_hat_type), top_hat_data.section_name)

	num_purlin_line_segments = length(purlin_line.inputs.segments)

	top_hat_cross_section_dimensions = Vector{NTuple{21, Float64}}(undef, num_purlin_line_segments)
	
	for i = 1:num_purlin_line_segments  #for now use the same TopHat everywhere, set up so that TopHat can change in the future
		
		top_hat_cross_section_dimensions[i] = tuple([top_hat_data[top_hat_section_index, :][i] for i=2:22]...)

	end

	# Define the TopHat material properties.
	top_hat_material_properties = [(29500.0, 0.30, 55.0, 70.0)]; #E, ν, Fy, Fu


	deck_index = findfirst(==(existing_deck_type), existing_deck_data.deck_name)
	
	# Define the TopHat punchout dimensions.  
	top_hat_punch_out_dimensions = (existing_deck_data[deck_index, 7], existing_deck_data[deck_index, 8]); #length, height

  	#Define the new deck details.

	new_deck_index = findfirst(==(new_deck_type), new_deck_data.deck_name)


	if !ismissing(new_deck_data[new_deck_index, 3]) #screw_fastened

		new_roof_panel_details = ("screw-fastened", new_deck_data[new_deck_index, 2], new_deck_data[new_deck_index, 3], new_deck_data[new_deck_index, 4], new_deck_data[new_deck_index, 5])

	elseif ismissing(new_deck_data[new_deck_index, 3]) #SSR
	
		new_roof_panel_details = ("vertical leg standing seam", new_deck_data[new_deck_index, 7], new_deck_data[new_deck_index, 6], 0.0, 0.0)

	end

	new_roof_panel_material_properties = (29500.0, 0.30, 55.0, 70.0)

	# Assemble the purlin line model, now with the addition of TopHat Framing.

	top_hat_purlin_line = TopHatDesigner.define(purlin_line.inputs.design_code, purlin_line.inputs.segments, purlin_line.inputs.spacing, purlin_line.inputs.roof_slope, purlin_line.inputs.cross_section_dimensions, top_hat_cross_section_dimensions, top_hat_punch_out_dimensions, purlin_line.inputs.material_properties, top_hat_material_properties, purlin_line.inputs.deck_details, purlin_line.inputs.deck_material_properties, new_roof_panel_details, new_roof_panel_material_properties, purlin_line.inputs.frame_flange_width, purlin_line.inputs.support_locations, purlin_line.inputs.purlin_frame_connections, purlin_line.inputs.bridging_locations)

	# Run a test to calculate the expected roof system failure pressure including the TopHat.

	#Run a gravity test.
	top_hat_purlin_line_gravity = deepcopy(top_hat_purlin_line)
	top_hat_purlin_line_gravity.loading_direction = "gravity"
	top_hat_purlin_line_gravity = TopHatDesigner.capacity(top_hat_purlin_line_gravity)

	#Run an uplift test.
	top_hat_purlin_line_uplift = deepcopy(top_hat_purlin_line)
	top_hat_purlin_line_uplift.loading_direction = "uplift"
	top_hat_purlin_line_uplift = TopHatDesigner.capacity(top_hat_purlin_line_uplift)

	return top_hat_purlin_line_gravity, top_hat_purlin_line_uplift

end

function generate_purlin_geometry(t, xcoords_center, ycoords_center, roof_slope)

	center_nodes = [xcoords_center ycoords_center zeros(Float64, length(xcoords_center))]

	center_nodes_rotated = LinesCurvesNodes.rotate_nodes(center_nodes, rotation_axis = "z", rotation_center = [0.0, 0.0, 0.0], θ=atan(roof_slope))
	
	cross_section = [[xcoords_center[i], ycoords_center[i]] for i in eachindex(xcoords_center)]
	
	unit_node_normals = CrossSection.Tools.calculate_cross_section_unit_node_normals(cross_section)
	
	outside = CrossSection.Tools.get_coords_along_node_normals(cross_section, unit_node_normals, t/2)
	X = [outside[i][1] for i in eachindex(outside)]
	Y = [outside[i][2] for i in eachindex(outside)]
	out_nodes = [X Y zeros(Float64, length(xcoords_center))]
	out_nodes_rotated = LinesCurvesNodes.rotate_nodes(out_nodes, rotation_axis = "z", rotation_center = [0.0, 0.0, 0.0], θ=atan(roof_slope))
	
	inside = CrossSection.Tools.get_coords_along_node_normals(cross_section, unit_node_normals, -t/2)
	X = [inside[i][1] for i in eachindex(inside)]
	Y = [inside[i][2] for i in eachindex(inside)]
	in_nodes = [X Y zeros(Float64, length(xcoords_center))]
	in_nodes_rotated = LinesCurvesNodes.rotate_nodes(in_nodes, rotation_axis = "z", rotation_center = [0.0, 0.0, 0.0], θ=atan(roof_slope))

	return center_nodes_rotated, out_nodes_rotated, in_nodes_rotated
	
end

function plot_purlin_geometry(t, xcoords_center, ycoords_center, roof_slope)

	center_nodes, out_nodes, in_nodes = generate_purlin_geometry(t, xcoords_center, ycoords_center, roof_slope)

	plot(center_nodes[:,1], center_nodes[:,2], aspect_ratio=:equal, linecolor = :grey, legend=false)

	plot!(out_nodes[:,1], out_nodes[:,2], aspect_ratio=:equal, linecolor = :grey, legend=false)

	plot!(in_nodes[:,1], in_nodes[:,2], aspect_ratio=:equal, linecolor = :grey, legend=false)
	
end

function plot_top_hat_purlin_geometry(t, xcoords_center_purlin, ycoords_center_purlin, roof_slope, xcoords_top_hat_purlin_line, ycoords_top_hat_purlin_line, xcoords_center_top_hat, ycoords_center_top_hat)

	center_nodes_purlin, out_nodes_purlin, in_nodes_purlin = generate_purlin_geometry(t, xcoords_center_purlin, ycoords_center_purlin, roof_slope)

	plot(center_nodes_purlin[:,1], center_nodes_purlin[:,2], aspect_ratio=:equal, linecolor = :grey, legend=false)

	plot!(out_nodes_purlin[:,1], out_nodes_purlin[:,2], aspect_ratio=:equal, linecolor = :grey, legend=false)

	plot!(in_nodes_purlin[:,1], in_nodes_purlin[:,2], aspect_ratio=:equal, linecolor = :grey, legend=false)

	top_hat_bottom_flange_centerline_node = 47
	Δx = xcoords_top_hat_purlin_line[top_hat_bottom_flange_centerline_node]
	Δy = ycoords_top_hat_purlin_line[top_hat_bottom_flange_centerline_node]
	Δz = 0.0
	
	center_nodes_top_hat, out_nodes_top_hat, in_nodes_top_hat = generate_top_hat_geometry(t, xcoords_center_top_hat, ycoords_center_top_hat, roof_slope, Δx, Δy, Δz)

	plot!(center_nodes_top_hat[:,1], center_nodes_top_hat[:,2], aspect_ratio=:equal, linecolor = :grey, legend=false)

	plot!(out_nodes_top_hat[:,1], out_nodes_top_hat[:,2], aspect_ratio=:equal, linecolor = :grey, legend=false)

	plot!(in_nodes_top_hat[:,1], in_nodes_top_hat[:,2], aspect_ratio=:equal, linecolor = :grey, legend=false)

end	


function plot_net_section_top_hat_purlin_geometry(t, xcoords_center_purlin, ycoords_center_purlin, xcoords_center_top_hat, ycoords_center_top_hat, roof_slope, xcoords_center_top_hat_purlin, ycoords_center_top_hat_purlin, net_top_hat_purlin_node_geometry)

	center_nodes_purlin, out_nodes_purlin, in_nodes_purlin = generate_purlin_geometry(t, xcoords_center_purlin, ycoords_center_purlin, roof_slope)

	plot(center_nodes_purlin[:,1], center_nodes_purlin[:,2], aspect_ratio=:equal, linecolor = :grey, legend=false)

	plot!(out_nodes_purlin[:,1], out_nodes_purlin[:,2], aspect_ratio=:equal, linecolor = :grey, legend=false)

	plot!(in_nodes_purlin[:,1], in_nodes_purlin[:,2], aspect_ratio=:equal, linecolor = :grey, legend=false)

	top_hat_bottom_flange_centerline_node = 47
	Δx = xcoords_center_top_hat_purlin[top_hat_bottom_flange_centerline_node]
	Δy = ycoords_center_top_hat_purlin[top_hat_bottom_flange_centerline_node]
	Δz = 0.0
	
	center_nodes_top_hat, out_nodes_top_hat, in_nodes_top_hat = generate_top_hat_geometry(t, xcoords_center_top_hat, ycoords_center_top_hat, roof_slope, Δx, Δy, Δz)

	num_top_hat_purlin_net_section_nodes = size(net_top_hat_purlin_node_geometry, 1)

	num_purlin_nodes = size(xcoords_center_purlin, 1)

	num_top_hat_net_section_nodes = num_top_hat_purlin_net_section_nodes - num_purlin_nodes

	plot_index = Int(num_top_hat_net_section_nodes/2)

	plot!(center_nodes_top_hat[1:plot_index,1], center_nodes_top_hat[1:plot_index,2], aspect_ratio=:equal, linecolor = :grey, legend=false)

	plot!(out_nodes_top_hat[1:plot_index,1], out_nodes_top_hat[1:plot_index,2], aspect_ratio=:equal, linecolor = :grey, legend=false)

	plot!(in_nodes_top_hat[1:plot_index,1], in_nodes_top_hat[1:plot_index,2], aspect_ratio=:equal, linecolor = :grey, legend=false)

	plot!(center_nodes_top_hat[end-plot_index+1:end,1], center_nodes_top_hat[end-plot_index+1:end,2], aspect_ratio=:equal, linecolor = :grey, legend=false)

	plot!(out_nodes_top_hat[end-plot_index+1:end,1], out_nodes_top_hat[end-plot_index+1:end,2], aspect_ratio=:equal, linecolor = :grey, legend=false)

	plot!(in_nodes_top_hat[end-plot_index+1:end,1], in_nodes_top_hat[end-plot_index+1:end,2], aspect_ratio=:equal, linecolor = :grey, legend=false)

end	



function generate_top_hat_geometry(t, xcoords_center, ycoords_center, roof_slope, Δx, Δy, Δz)

	center_nodes = [xcoords_center ycoords_center zeros(Float64, length(xcoords_center))]

	center_nodes = LinesCurvesNodes.shift_nodes(center_nodes, Δx=Δx, Δy=Δy, Δz=Δz)

	center_nodes_rotated = LinesCurvesNodes.rotate_nodes(center_nodes, rotation_axis = "z", rotation_center = [0.0, 0.0, 0.0], θ=atan(roof_slope))
	

	cross_section = [[xcoords_center[i], ycoords_center[i]] for i in eachindex(xcoords_center)]
	
	unit_node_normals = CrossSection.Tools.calculate_cross_section_unit_node_normals(cross_section)
	
	outside = CrossSection.Tools.get_coords_along_node_normals(cross_section, unit_node_normals, t/2)
	X = [outside[i][1] for i in eachindex(outside)]
	Y = [outside[i][2] for i in eachindex(outside)]
	out_nodes = [X Y zeros(Float64, length(xcoords_center))]
	out_nodes = LinesCurvesNodes.shift_nodes(out_nodes, Δx=Δx, Δy=Δy, Δz=Δz)
	out_nodes_rotated = LinesCurvesNodes.rotate_nodes(out_nodes, rotation_axis = "z", rotation_center = [0.0, 0.0, 0.0], θ=atan(roof_slope))
	
	inside = CrossSection.Tools.get_coords_along_node_normals(cross_section, unit_node_normals, -t/2)
	X = [inside[i][1] for i in eachindex(inside)]
	Y = [inside[i][2] for i in eachindex(inside)]
	in_nodes = [X Y zeros(Float64, length(xcoords_center))]
	in_nodes = LinesCurvesNodes.shift_nodes(in_nodes, Δx=Δx, Δy=Δy, Δz=Δz)
	in_nodes_rotated = LinesCurvesNodes.rotate_nodes(in_nodes, rotation_axis = "z", rotation_center = [0.0, 0.0, 0.0], θ=atan(roof_slope))


   ####3

	# center_nodes_rotated = Geometry.rotate_nodes(center_nodes, rotation_axis = "z", rotation_center = [0.0, 0.0, 0.0], θ=atan(roof_slope))
	
	# 	#Calculate surface normals.
	# closed_or_open = 1
	# unitnormals = SectionProperties.surface_normals(center_nodes[:, 1], center_nodes[:, 2], closed_or_open)

	
	# nodenormals = SectionProperties.avg_node_normals(unitnormals, closed_or_open)


	# xcoords_out, ycoords_out = SectionProperties.xycoords_along_normal(center_nodes[:, 1], center_nodes[:, 2], nodenormals, t/2)

	# out_nodes_rotated = [xcoords_out ycoords_out zeros(Float64, length(xcoords_center))]

	# # out_nodes_rotated = Geometry.rotate_nodes(out_nodes, rotation_axis = "z", rotation_center = [0.0, 0.0, 0.0], θ=atan(roof_slope))

	# xcoords_in, ycoords_in = SectionProperties.xycoords_along_normal(center_nodes[:, 1], center_nodes[:, 2], nodenormals, -t/2)

	# in_nodes_rotated = [xcoords_in ycoords_in zeros(Float64, length(xcoords_center))]

	# in_nodes_rotated = Geometry.rotate_nodes(in_nodes, rotation_axis = "z", rotation_center = [0.0, 0.0, 0.0], θ=atan(roof_slope))

	return center_nodes_rotated, out_nodes_rotated, in_nodes_rotated
	
end

end
