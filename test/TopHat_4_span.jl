using TopHatDesigner, CSV, DataFrames


purlin_data = CSV.read("database/Purlins.csv",
DataFrame);

top_hat_data = CSV.read("database/TopHats.csv",
DataFrame);

existing_deck_data = CSV.read("database/Existing_Deck.csv",
DataFrame);

new_deck_data = CSV.read("database/New_Deck.csv",
DataFrame);

purlin_type_1 = "Z8x2.5 067"

purlin_type_2 = "none"

purlin_spans = (26.0, 27.0, 27.0, 26.0)  #ft

purlin_size_span_assignment = (1, 1, 1, 1)

purlin_laps = (2.0, 2.0, 2.0, 2.0, 2.0, 2.0)

purlin_spacing = 4.0  #ft

frame_flange_width = 16.0  #in

roof_slope = 1/12

existing_deck_type = "SSR Ultra-Dek 24 in. 24 ga "

purlin_line, purlin_line_uplift = UI.existing_roof_UI_mapper(purlin_spans, purlin_laps, purlin_spacing, roof_slope, purlin_data, existing_deck_type, existing_deck_data, frame_flange_width, purlin_type_1, purlin_type_2, purlin_size_span_assignment);

purlin_line.applied_pressure*1000*144

purlin_line.failure_location

purlin_line.failure_limit_state

purlin_line_uplift.applied_pressure*1000*144

using Plots
plot(purlin_line.model.z, purlin_line.internal_forces.T)

top_hat_type = "TH3.5 071"
		
new_deck_type = "SSR SuperLok-Dek 16 in. 24 ga "


top_hat_purlin_line_gravity, top_hat_purlin_line_uplift = TopHatDesigner.UI.retrofit_UI_mapper(purlin_line, top_hat_data, top_hat_type, existing_deck_type, existing_deck_data, new_deck_type, new_deck_data)

top_hat_purlin_line_gravity.applied_pressure*1000*144

top_hat_purlin_line_uplift.applied_pressure*1000*144

purlin_line.failure_limit_state

top_hat_purlin_line_gravity.Β_distortional_gradient_factor

using Plots
plot(top_hat_purlin_line_gravity.model.z, top_hat_purlin_line_gravity.biaxial_bending_demand_to_capacity.action_Mxx)
plot!(top_hat_purlin_line_gravity.model.z, top_hat_purlin_line_gravity.biaxial_bending_demand_to_capacity.action_Myy)
plot!(top_hat_purlin_line_gravity.model.z, top_hat_purlin_line_gravity.biaxial_bending_demand_to_capacity.action_Myy)

