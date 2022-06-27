using TopHatDesigner, CSV, DataFrames, LinesCurvesNodes, CrossSection


purlin_data = CSV.read("/Users/crismoen/.julia/dev/TopHatDesigner/database/Purlins.csv",
DataFrame);

top_hat_data = CSV.read("/Users/crismoen/.julia/dev/TopHatDesigner/database/TopHats.csv",
DataFrame);

existing_deck_data = CSV.read("/Users/crismoen/.julia/dev/TopHatDesigner/database/Existing_Deck.csv",
DataFrame);

new_deck_data = CSV.read("/Users/crismoen/.julia/dev/TopHatDesigner/database/New_Deck.csv",
DataFrame);

purlin_type_1 = "C8x2.5 060"

purlin_type_2 = "none"

purlin_spans = (25.0)  #ft

purlin_size_span_assignment = (1)

purlin_laps = ()

purlin_spacing = 4.0  #ft

frame_flange_width = 16.0  #in

roof_slope = 1/12

existing_deck_type = "SSR Ultra-Dek 24 in. 24 ga "

purlin_line, purlin_line_uplift = UI.existing_roof_UI_mapper(purlin_spans, purlin_laps, purlin_spacing, roof_slope, purlin_data, existing_deck_type, existing_deck_data, frame_flange_width, purlin_type_1, purlin_type_2, purlin_size_span_assignment);


top_hat_type = "TH3.5 071"
		
new_deck_type = "SSR SuperLok-Dek 16 in. 24 ga "

top_hat_purlin_line, top_hat_purlin_line_uplift = TopHatDesigner.UI.retrofit_UI_mapper(purlin_line, top_hat_data, top_hat_type, existing_deck_type, existing_deck_data, new_deck_type, new_deck_data);


UI.plot_purlin_geometry(purlin_line.inputs.cross_section_dimensions[1][2], purlin_line.cross_section_data[1].node_geometry[:,1], purlin_line.cross_section_data[1].node_geometry[:,2], roof_slope)


UI.plot_top_hat_purlin_geometry(top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[1][1], top_hat_purlin_line.purlin_cross_section_data[1].node_geometry[:,1], top_hat_purlin_line.purlin_cross_section_data[1].node_geometry[:,2], roof_slope, top_hat_purlin_line.top_hat_purlin_cross_section_data[1].node_geometry[:,1], top_hat_purlin_line.top_hat_purlin_cross_section_data[1].node_geometry[:,2], top_hat_purlin_line.top_hat_cross_section_data[1].node_geometry[:,1], top_hat_purlin_line.top_hat_cross_section_data[1].node_geometry[:,2])

UI.plot_net_section_top_hat_purlin_geometry(purlin_line.inputs.cross_section_dimensions[1][2], top_hat_purlin_line.purlin_cross_section_data[1].node_geometry[:,1], top_hat_purlin_line.purlin_cross_section_data[1].node_geometry[:,2], top_hat_purlin_line.top_hat_cross_section_data[1].node_geometry[:,1], top_hat_purlin_line.top_hat_cross_section_data[1].node_geometry[:,2], roof_slope, top_hat_purlin_line.top_hat_purlin_cross_section_data[1].node_geometry[:,1],  top_hat_purlin_line.top_hat_purlin_cross_section_data[1].node_geometry[:,2], top_hat_purlin_line.top_hat_purlin_net_cross_section_data[1].node_geometry)




