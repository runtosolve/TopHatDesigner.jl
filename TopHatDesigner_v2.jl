### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ e37326d5-d9c8-4512-960f-82d912c00273
begin
    import Pkg
    Pkg.activate()

    using PlutoUI, Images, Dates, CSV, DataFrames, PurlinLine, TopHatDesigner, StructuresKit
	
end

# ╔═╡ 43ba7d07-7cd6-4803-ae5c-6b23a0cf023b
purlin_data = CSV.read("/Volumes/GoogleDrive/Shared drives/RunToSolve/Projects/TopHat Framing/TopHatDesigner_UI/database/Purlins.csv",
                             DataFrame);

# ╔═╡ 91637826-04f4-488e-9ee5-024dc5221cdc
top_hat_data = CSV.read("/Volumes/GoogleDrive/Shared drives/RunToSolve/Projects/TopHat Framing/TopHatDesigner_UI/database/TopHats.csv",
                             DataFrame);

# ╔═╡ 299130e4-be9a-4b5c-a851-8bf3aa37d28b
existing_deck_data = CSV.read("/Volumes/GoogleDrive/Shared drives/RunToSolve/Projects/TopHat Framing/TopHatDesigner_UI/database/Existing_Deck.csv",
                             DataFrame);

# ╔═╡ 7c9236f3-7a63-4686-bd87-cde4efb58e87
new_deck_data = CSV.read("/Volumes/GoogleDrive/Shared drives/RunToSolve/Projects/TopHat Framing/TopHatDesigner_UI/database/New_Deck.csv",
                             DataFrame);

# ╔═╡ 2a2939b5-b0de-4b51-9a14-3fdf5d439b93
load("/Volumes/GoogleDrive/Shared drives/RunToSolve/Projects/TopHat Framing/TopHatDesigner_UI/tophat-horz.png")

# ╔═╡ c2353bb3-ee8c-4d55-9447-470427c22b06
@bind project_details TextField((30,5); default="Project details")

# ╔═╡ 96f90537-0b4d-4d48-927b-01492e3789ef
@bind report_date DateField(default=today())

# ╔═╡ 99299f0c-30ee-4807-a7a2-d4509b4680ab
md" ## Build existing roof system."

# ╔═╡ 5d180e53-27aa-4bb2-9ab1-81cc2737ab3b
purlin_spans = (25.0)  #ft

# ╔═╡ 1d9b00aa-7f6b-4f7e-9da8-e1f0b1ace647
md"Purlin size $(@bind purlin_type Select(purlin_data[:, 1]))"

# ╔═╡ 710a97bb-6cd1-457c-a352-23428408de55
purlin_laps = ()

# ╔═╡ dde7f4c2-2212-4244-a78b-8fe12b6c8d0e
purlin_spacing = 4.0  #ft

# ╔═╡ 1a1727db-c828-4317-8ccd-2be491ee48c0
frame_flange_width = 16.0  #in

# ╔═╡ 9eebfba0-7913-40fd-bda5-b1d1178c741a
roof_slope = 1/12

# ╔═╡ e20e6735-ae8a-4ae0-99bb-f563a602afbc
md" Existing roof deck type $(@bind existing_deck_type Select(existing_deck_data[:, 1]))"

# ╔═╡ a2083ad0-b14d-4968-9821-207b4af24121
	section_index = findfirst(==(purlin_type), purlin_data.section_name)

# ╔═╡ dd6bbcfe-b781-4bf5-8485-1c4b25380ebc
md" ## Calculate existing roof system strength."

# ╔═╡ bc141195-cc2e-437d-907c-92ba11e4d55f
function existing_roof_UI_mapper(purlin_spans, purlin_laps, purlin_spacing, roof_slope, purlin_type, purlin_data, existing_deck_type, existing_deck_data, frame_flange_width)

	design_code = "AISI S100-16 ASD";

	num_purlin_segments = length(purlin_spans) + length(purlin_laps)	

	purlin_segments = Array{Tuple{Float64, Int64, Int64}, 1}(undef, num_purlin_segments)

	for i in eachindex(purlin_spans)

		purlin_segments[i] = (purlin_spans[i]*12, 1, 1)

	end

	purlin_spacing = purlin_spacing * 12.0

	roof_slope = rad2deg(atan(roof_slope))

	section_index = findfirst(==(purlin_type), purlin_data.section_name)
	
	purlin_cross_section_dimensions = [tuple([purlin_data[section_index, :][i] for i=2:17]...)]

	purlin_material_properties = [(29500.0, 0.30, 50.0, 70.0)];  #E, ν, Fy, Fu

	deck_index = findfirst(==(existing_deck_type), existing_deck_data.deck_name)
	existing_roof_panel_details = ("screw-fastened", existing_deck_data[deck_index, 2], existing_deck_data[deck_index, 3], existing_deck_data[deck_index, 4], existing_deck_data[deck_index, 5])

	existing_roof_panel_material_properties = (29500.0, 0.30, 55.0, 70.0);  #E, ν, Fy, Fu

	support_locations = [0.0, purlin_spans[end]*12.0]

	purlin_frame_connections = "anti-roll clip"

	intermediate_bridging_locations = [ ]

		purlin_line = PurlinLine.build(design_code, purlin_segments, purlin_spacing, roof_slope, purlin_cross_section_dimensions, purlin_material_properties, existing_roof_panel_details, existing_roof_panel_material_properties, frame_flange_width, support_locations, purlin_frame_connections, intermediate_bridging_locations)

	#Run a gravity test.
	purlin_line.loading_direction = "gravity"
	purlin_line = PurlinLine.test(purlin_line)

	return purlin_line

end;

# ╔═╡ 155b41ff-d261-40df-807d-37027d220187
purlin_line = existing_roof_UI_mapper(purlin_spans, purlin_laps, purlin_spacing, roof_slope, purlin_type, purlin_data, existing_deck_type, existing_deck_data, frame_flange_width);

# ╔═╡ 40af7c07-ff81-4617-a220-9db2435247f2
begin

	xlims=(-9.0, 9.0);
	ylims = (0.0, 18.0);
	markershape = :none;
	StructuresKit.Visualize.show_multi_branch_cross_section(purlin_line.cross_section_data[1].node_geometry[:,1], purlin_line.cross_section_data[1].node_geometry[:,2], purlin_line.cross_section_data[1].element_definitions, markershape, xlims, ylims)

end

# ╔═╡ d2246077-5cf8-4f00-81af-9e5922bec619
md"**Existing roof system strength = $(round(purlin_line.applied_pressure*1000*144, digits=1)) psf**"

# ╔═╡ f67cf06d-2fe6-401f-ab71-1ecea9aa373f
md" ## Add TopHat and new roof."

# ╔═╡ 6cebdfda-7d1c-4df6-9e64-b6afad1c7f6d
md" TopHat size $(@bind top_hat_type Select(top_hat_data[:, 1]))"
		

# ╔═╡ bfac357a-4b13-471e-84e8-edf4272065ba
md" New roof deck type $(@bind new_deck_type Select(new_deck_data[:, 1]))"

# ╔═╡ e8b4678d-a329-4403-bb9e-00f827a86c98
new_deck_index = findfirst(==(new_deck_type), new_deck_data.deck_name)

# ╔═╡ b4576ffe-199e-4763-a542-a015fe1d4164
new_roof_panel_details = ("vertical leg standing seam", new_deck_data[new_deck_index, 6], new_deck_data[new_deck_index, 7], 0.0, 0.0)

# ╔═╡ 9ce26903-de7d-49b8-a7f8-7e98c892451a
ismissing(new_deck_data[new_deck_index, 3])

# ╔═╡ 4d47c72c-52a9-4caa-9516-446296e73443
md" ## Calculate retrofitted roof system strength."

# ╔═╡ 2df480a5-d4c2-4d8a-a78e-4952b23939f1
function retrofit_UI_mapper(purlin_line, top_hat_data, top_hat_type, existing_deck_type, existing_deck_data, new_deck_type, new_deck_data)


	top_hat_section_index = findfirst(==(top_hat_type), top_hat_data.section_name)
	
	top_hat_cross_section_dimensions = [tuple([top_hat_data[top_hat_section_index, :][i] for i=2:22]...)]

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
	top_hat_purlin_line.loading_direction = "gravity";
	top_hat_purlin_line = TopHatDesigner.capacity(top_hat_purlin_line)

	return top_hat_purlin_line

end;

# ╔═╡ 33f14b20-ecce-43e9-884b-629ef5a9ac40
top_hat_purlin_line = retrofit_UI_mapper(purlin_line, top_hat_data, top_hat_type, existing_deck_type, existing_deck_data, new_deck_type, new_deck_data);

# ╔═╡ 828e6235-3b5f-47f0-a004-62b00bf14ff9
begin
	StructuresKit.Visualize.show_multi_branch_cross_section(top_hat_purlin_line.top_hat_purlin_cross_section_data[1].node_geometry[:,1], top_hat_purlin_line.top_hat_purlin_cross_section_data[1].node_geometry[:,2], top_hat_purlin_line.top_hat_purlin_cross_section_data[1].element_definitions, markershape, xlims, ylims)

end

# ╔═╡ 8fe61b54-e9e1-4aa6-be50-b0812d86a1cf
md"**Retrofitted roof system strength = $(round(top_hat_purlin_line.applied_pressure*1000*144, digits=1)) psf**"

# ╔═╡ Cell order:
# ╠═e37326d5-d9c8-4512-960f-82d912c00273
# ╠═43ba7d07-7cd6-4803-ae5c-6b23a0cf023b
# ╠═91637826-04f4-488e-9ee5-024dc5221cdc
# ╠═299130e4-be9a-4b5c-a851-8bf3aa37d28b
# ╠═7c9236f3-7a63-4686-bd87-cde4efb58e87
# ╟─2a2939b5-b0de-4b51-9a14-3fdf5d439b93
# ╟─c2353bb3-ee8c-4d55-9447-470427c22b06
# ╟─96f90537-0b4d-4d48-927b-01492e3789ef
# ╟─99299f0c-30ee-4807-a7a2-d4509b4680ab
# ╠═5d180e53-27aa-4bb2-9ab1-81cc2737ab3b
# ╠═1d9b00aa-7f6b-4f7e-9da8-e1f0b1ace647
# ╠═710a97bb-6cd1-457c-a352-23428408de55
# ╠═dde7f4c2-2212-4244-a78b-8fe12b6c8d0e
# ╠═1a1727db-c828-4317-8ccd-2be491ee48c0
# ╠═9eebfba0-7913-40fd-bda5-b1d1178c741a
# ╠═e20e6735-ae8a-4ae0-99bb-f563a602afbc
# ╠═a2083ad0-b14d-4968-9821-207b4af24121
# ╟─dd6bbcfe-b781-4bf5-8485-1c4b25380ebc
# ╟─bc141195-cc2e-437d-907c-92ba11e4d55f
# ╟─155b41ff-d261-40df-807d-37027d220187
# ╠═40af7c07-ff81-4617-a220-9db2435247f2
# ╟─d2246077-5cf8-4f00-81af-9e5922bec619
# ╟─f67cf06d-2fe6-401f-ab71-1ecea9aa373f
# ╠═6cebdfda-7d1c-4df6-9e64-b6afad1c7f6d
# ╠═bfac357a-4b13-471e-84e8-edf4272065ba
# ╠═e8b4678d-a329-4403-bb9e-00f827a86c98
# ╠═b4576ffe-199e-4763-a542-a015fe1d4164
# ╠═9ce26903-de7d-49b8-a7f8-7e98c892451a
# ╟─4d47c72c-52a9-4caa-9516-446296e73443
# ╠═2df480a5-d4c2-4d8a-a78e-4952b23939f1
# ╟─33f14b20-ecce-43e9-884b-629ef5a9ac40
# ╟─828e6235-3b5f-47f0-a004-62b00bf14ff9
# ╟─8fe61b54-e9e1-4aa6-be50-b0812d86a1cf
