module TopHatDesigner

using PurlinLine
using StructuresKit

using NumericalIntegration

export define, analysis, capacity

struct Inputs

    design_code::String
    segments::Vector{Tuple{Float64, Int64, Int64}}
    spacing::Float64
    roof_slope::Float64
    purlin_cross_section_dimensions::Vector{Tuple{String, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64}}
    top_hat_cross_section_dimensions::Vector{Tuple{Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64}}
    top_hat_punch_out_dimensions::NTuple{2, Float64}
    purlin_material_properties::Vector{NTuple{4, Float64}}
    top_hat_material_properties::Vector{NTuple{4, Float64}}
    new_deck_details::Tuple{String, Float64, Float64, Float64, Float64}
    new_deck_material_properties::NTuple{4, Float64}
    frame_flange_width::Float64
    support_locations::Vector{Float64}
    purlin_frame_connections::String
    bridging_locations::Vector{Float64}

end

mutable struct TopHatDesignerObject

    inputs::TopHatDesigner.Inputs

    applied_pressure::Float64

    loading_direction::String

    purlin_cross_section_data::Array{PurlinLine.CrossSectionData}

    free_flange_cross_section_data::Array{PurlinLine.CrossSectionData}

    top_hat_cross_section_data::Array{PurlinLine.CrossSectionData}

    top_hat_purlin_cross_section_data::Array{PurlinLine.CrossSectionData}

    top_hat_purlin_net_cross_section_data::Array{PurlinLine.CrossSectionData}

    top_hat_purlin_distortional_net_cross_section_data::Array{PurlinLine.CrossSectionData}

    bracing_data::Array{PurlinLine.BracingData}

    new_deck_bracing_data::Array{PurlinLine.BracingData}

    free_flange_data::Array{PurlinLine.FreeFlangeData}

    local_buckling_xx_pos::Array{PurlinLine.ElasticBucklingData}
    local_buckling_xx_net_pos::Array{PurlinLine.ElasticBucklingData}
    local_buckling_xx_neg::Array{PurlinLine.ElasticBucklingData}
    local_buckling_yy_pos::Array{PurlinLine.ElasticBucklingData}
    local_buckling_yy_neg::Array{PurlinLine.ElasticBucklingData}
    distortional_buckling_xx_pos::Array{PurlinLine.ElasticBucklingData}
    distortional_buckling_xx_net_pos::Array{PurlinLine.ElasticBucklingData}
    distortional_buckling_xx_neg::Array{PurlinLine.ElasticBucklingData}

    yielding_flexural_strength_xx::Array{PurlinLine.YieldingFlexuralStrengthData}
    yielding_flexural_strength_xx_net::Array{PurlinLine.YieldingFlexuralStrengthData}
    yielding_flexural_strength_yy::Array{PurlinLine.YieldingFlexuralStrengthData}
    yielding_flexural_strength_free_flange_yy::Array{PurlinLine.YieldingFlexuralStrengthData}

    local_global_flexural_strength_xx_no_hole::Array{PurlinLine.LocalGlobalFlexuralStrengthData}
    local_global_flexural_strength_xx_hole::Array{PurlinLine.LocalGlobalFlexuralStrengthData}
    local_global_flexural_strength_xx::Array{PurlinLine.LocalGlobalFlexuralStrengthData}
    local_global_flexural_strength_yy::Array{PurlinLine.LocalGlobalFlexuralStrengthData}
    local_global_flexural_strength_free_flange_yy::Array{PurlinLine.LocalGlobalFlexuralStrengthData}

    distortional_flexural_strength_xx::Array{PurlinLine.DistortionalFlexuralStrengthData}

    torsion_strength::Array{PurlinLine.TorsionStrengthData}

    shear_strength_purlin::Array{PurlinLine.ShearStrengthData}
    shear_strength_top_hat::Array{PurlinLine.ShearStrengthData}
    shear_strength::Array{PurlinLine.ShearStrengthData}

    web_crippling::Array{PurlinLine.WebCripplingData}

    model::ThinWalledBeam.Model

    free_flange_model::BeamColumn.Model

    internal_forces::PurlinLine.InternalForceData

    free_flange_internal_forces::PurlinLine.InternalForceData

    support_reactions::PurlinLine.Reactions

    flexure_torsion_demand_to_capacity::PurlinLine.FlexureTorsion_DemandToCapacity_Data
    biaxial_bending_demand_to_capacity::PurlinLine.BiaxialBending_DemandToCapacity_Data
    distortional_demand_to_capacity::Array{Float64}
    flexure_shear_demand_to_capacity::Array{Float64}
    web_crippling_demand_to_capacity::Array{Float64}

    expected_strengths::PurlinLine.ExpectedStrengths

    failure_limit_state::String

    failure_location::Float64

    TopHatDesignerObject() = new()

end


function define_top_hat_cross_sections(cross_section_dimensions, n, n_radius)

    num_top_hat_sections = size(cross_section_dimensions)[1]

    cross_section_data = Vector{PurlinLine.CrossSectionData}(undef, num_top_hat_sections)

    for i = 1:num_top_hat_sections

        #Map dimensions to cross-section nomenclature.
        t = cross_section_dimensions[i][1]
        d_top_minus_y = cross_section_dimensions[i][2]
        b_top_minus_y = cross_section_dimensions[i][3]
        h_minus_y = cross_section_dimensions[i][4]
        b_bottom = cross_section_dimensions[i][5]
        h_plus_y = cross_section_dimensions[i][6]
        b_top_plus_y = cross_section_dimensions[i][7]
        d_top_plus_y = cross_section_dimensions[i][8]
        α1 = cross_section_dimensions[i][9]
        α2 = cross_section_dimensions[i][10]
        α3 = cross_section_dimensions[i][11]
        α4 = cross_section_dimensions[i][12]
        α5 = cross_section_dimensions[i][13]
        α6 = cross_section_dimensions[i][14]
        α7 = cross_section_dimensions[i][15]
        r1 = cross_section_dimensions[i][16]
        r2 = cross_section_dimensions[i][17]
        r3 = cross_section_dimensions[i][18]
        r4 = cross_section_dimensions[i][19]
        r5 = cross_section_dimensions[i][20]
        r6 = cross_section_dimensions[i][21]

        #Define straight-line lengths on the top cross-section surface.   
        ΔL = [d_top_minus_y, b_top_minus_y, h_minus_y - t, b_bottom - 2*t, h_plus_y - t, b_top_plus_y, d_top_plus_y]
        θ = [α1, α2, α3, α4, α5, α6, α7]

        #Note that the outside radius is used at the top flanges, and the inside radius is used for the bottom flange.
        radius = [r1, r2, r3-t, r4-t, r5, r6]

        closed_or_open = 1

        top_hat_section = CrossSection.Feature(ΔL, θ, n, radius, n_radius, closed_or_open)

        #Calculate the out-to-out surface coordinates.
        xcoords_out, ycoords_out = CrossSection.get_xy_coordinates(top_hat_section)

        #Calculate centerline coordinates.
        unitnormals = CrossSection.surface_normals(xcoords_out, ycoords_out, closed_or_open)
        nodenormals = CrossSection.avg_node_normals(unitnormals, closed_or_open)
        xcoords_center, ycoords_center = CrossSection.xycoords_along_normal(xcoords_out, ycoords_out, nodenormals, -t/2)

        #Shift y coordinates so that the bottom face is at y = 0.
        ycoords_center = ycoords_center .- minimum(ycoords_center) .+ t/2

        #Shift x coordinates so that the section centerline is at x = 0.
        index = n[1] + n[2] + n[3] + floor(Int, n[4]/2) + n_radius[1] + n_radius[2] + n_radius[3] + 1
        xcoords_center = xcoords_center .- xcoords_center[index]

        #Package nodal geometry.
        node_geometry = [xcoords_center ycoords_center]

        #Define cross-section element connectivity and thicknesses.
        num_cross_section_nodes = length(xcoords_center)
        element_info = [1:(num_cross_section_nodes - 1) 2:num_cross_section_nodes ones(num_cross_section_nodes - 1) * t]

        #Calculate section properties.
        section_properties = CUFSM.cutwp_prop2(node_geometry, element_info)

        #Add cross section information to data structure.
        cross_section_data[i] = PurlinLine.CrossSectionData(n, n_radius, node_geometry, element_info, section_properties, nothing)

    end

    return cross_section_data

end

function combine_top_hat_purlin_geometry(purlin_cross_section_dimensions, purlin_cross_section_data, top_hat_cross_section_data)

    #Find purlin top flange centerline.
    purlin_top_flange_centerline_index = sum(purlin_cross_section_data.n[1:3]) + floor(Int,purlin_cross_section_data.n[4]/2) + sum(purlin_cross_section_data.n_radius[1:3]) + 1
    purlin_top_flange_centerline_geometry = purlin_cross_section_data.node_geometry[purlin_top_flange_centerline_index, :]

    #Shift TopHat geometry to sit on top of purlin top flange.
    top_hat_node_geometry = deepcopy(top_hat_cross_section_data.node_geometry)
    purlin_t = purlin_cross_section_dimensions[2]

    top_hat_node_geometry[:, 1] = top_hat_node_geometry[:, 1] .+ purlin_top_flange_centerline_geometry[1]
    top_hat_node_geometry[:, 2] = top_hat_node_geometry[:, 2] .+ purlin_top_flange_centerline_geometry[2] .+ purlin_t/2 

    # #Shift TopHat node numbers after purlin node numbers.

    num_purlin_nodes = size(purlin_cross_section_data.node_geometry)[1]
    top_hat_element_definitions = deepcopy(top_hat_cross_section_data.element_definitions)

    top_hat_element_definitions[:,1:2] = top_hat_element_definitions[:,1:2] .+ num_purlin_nodes

    #Combine node geometry and element definitions for TopHat and purlin.

    top_hat_purlin_node_geometry = [purlin_cross_section_data.node_geometry; top_hat_node_geometry]

    top_hat_purlin_element_definitions = [purlin_cross_section_data.element_definitions; top_hat_element_definitions]


    return top_hat_purlin_node_geometry, top_hat_purlin_element_definitions

end



function define_top_hat_purlin_cross_sections(purlin_cross_section_dimensions, purlin_cross_section_data, top_hat_cross_section_data, purlin_plastic_cross_section_data, top_hat_plastic_cross_section_data)

    #Assume that for each purlin segment, there is a purlin cross-section and a TopHat cross-section defined.
    num_strengthened_sections = size(purlin_cross_section_data)[1]

    top_hat_purlin_cross_section_data = Vector{PurlinLine.CrossSectionData}(undef, num_strengthened_sections)

    for i=1:num_strengthened_sections

        top_hat_purlin_node_geometry, top_hat_purlin_element_definitions = combine_top_hat_purlin_geometry(purlin_cross_section_dimensions[i], purlin_cross_section_data[i], top_hat_cross_section_data[i])

        #Calculate elastic section properties.
        top_hat_purlin_section_properties = CUFSM.cutwp_prop2(top_hat_purlin_node_geometry, top_hat_purlin_element_definitions)

        #Combine discretization info.
        n = [purlin_cross_section_data[i].n; top_hat_cross_section_data[i].n]
        n_radius = [purlin_cross_section_data[i].n_radius; top_hat_cross_section_data[i].n_radius]

        #######
        #Calculate TopHat+purlin plastic neutral axis and plastic modulus.
        top_hat_purlin_plastic_node_geometry, top_hat_purlin_plastic_element_definitions = combine_top_hat_purlin_geometry(purlin_cross_section_dimensions[i], purlin_plastic_cross_section_data[i], top_hat_plastic_cross_section_data[i])

        about_axis = "x"  #The strong axis plastic properties are needed for now.  
        top_hat_purlin_plastic_section_properties = StructuresKit.CrossSection.calculate_plastic_section_properties(top_hat_purlin_plastic_node_geometry, top_hat_purlin_plastic_element_definitions, about_axis)

        #Add cross section information to data structure.
        top_hat_purlin_cross_section_data[i] = PurlinLine.CrossSectionData(n, n_radius, top_hat_purlin_node_geometry, top_hat_purlin_element_definitions, top_hat_purlin_section_properties, top_hat_purlin_plastic_section_properties)

    end

    return top_hat_purlin_cross_section_data

end

function define_new_deck_bracing_properties(top_hat_purlin_line)

    num_purlin_segments = size(top_hat_purlin_line.inputs.segments)[1]

    bracing_data = Array{PurlinLine.BracingData, 1}(undef, num_purlin_segments)

    if top_hat_purlin_line.inputs.new_deck_details[1] == "screw-fastened"

        #Define the deck to TopHat screw-fastened connection spacing.
        deck_top_hat_fastener_spacing = top_hat_purlin_line.inputs.new_deck_details[3]

        #Define the deck to purlin screw diameter.
        deck_top_hat_fastener_diameter = top_hat_purlin_line.inputs.new_deck_details[4]

        #Define the nominal shear strength of the typical screw.
        Fss = top_hat_purlin_line.inputs.new_deck_details[5]

        #Define the roof deck base metal thickness.
        t_roof_deck = top_hat_purlin_line.inputs.new_deck_details[2]

        #Define roof deck steel elastic modulus.
        E_roof_deck = top_hat_purlin_line.inputs.new_deck_material_properties[1]

        #Define roof deck steel ultimate yield stress.
        Fu_roof_deck = top_hat_purlin_line.inputs.new_deck_material_properties[4]

        #Define the distance between fasteners as the distortional discrete bracing length.
        Lm = deck_top_hat_fastener_spacing

        #Loop over all the purlin segments in the line.  
        #Assume that there is a TopHat segment for every purlin segment.
        for i = 1:num_purlin_segments

            #Define the section property index associated with purlin segment i.
            section_index = top_hat_purlin_line.inputs.segments[i][2]

            #Define the material property index associated with purlin segment i.
            material_index = top_hat_purlin_line.inputs.segments[i][3]

            #Define TopHat steel elastic modulus.
            E_top_hat = top_hat_purlin_line.inputs.top_hat_material_properties[material_index][1]

            #Define TopHat steel Poisson's ratio.
            μ_top_hat = top_hat_purlin_line.inputs.top_hat_material_properties[material_index][2]

            #Define TopHat steel ultimate stress.
            Fu_top_hat = top_hat_purlin_line.inputs.top_hat_material_properties[material_index][4]

            #Define the TopHat top flange width.
            b_top = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[section_index][3]

            #Define TopHat base metal thickness.
            t_top_hat = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[section_index][1]

            #Define out-to-out TopHat web depth.
            ho = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[section_index][4]

            #Define TopHat top flange lip length.
            d_top = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[section_index][2]

            #Define TopHat top flange lip angle from the horizon, in degrees.
            θ_top = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[section_index][9]

            #Define the location from the TopHat top flange pivot point to the fastener.  Assume the fastener is centered in the flange.
            c = b_top/2

            #Define the deck fastener pull-through plate stiffness.  Assume the fastener is centered between two panel ribs.
            kp = PurlinLine.deck_pull_through_fastener_stiffness(top_hat_purlin_line.inputs.new_deck_material_properties, b_top, t_roof_deck)

            #Apply Cee or Zee binary.   Assume the TopHat behaves like a Z for this stiffness calculation.
            CorZ = 1
      
            #Calculate the rotational stiffness provided to each TopHat flange by the screw-fastened connection between the deck and the TopHat.  It is assumed that the deck flexural stiffness is much higher than the connection stiffness.
            kϕ = Connections.cfs_rot_screwfastened_k(b_top, c, deck_top_hat_fastener_spacing, t_top_hat, kp, E_top_hat, CorZ)

            #Calculate the TopHat distortional buckling half-wavelength.

            #Calculate top flange + lip section properties.
            Af, Jf, Ixf, Iyf, Ixyf, Cwf, xof,  hxf, hyf, yof = AISIS10016.table23131(CorZ, t_top_hat, b_top, d_top, θ_top)

            #Calculate the TopHat distortional buckling half-wavelength.
            Lcrd, L = AISIS10016.app23334(ho, μ_top_hat, t_top_hat, Ixf, xof, hxf, Cwf, Ixyf, Iyf, Lm)

            #If Lcrd is longer than the fastener spacing, then the distortional buckling will be restrained by the deck.
            if Lcrd >= Lm
                kϕ_dist = kϕ
            else
                kϕ_dist = 0.0
            end

            #Approximate the lateral stiffness provided to the top of each TopHat flange by the screw-fastened connection between the deck and the TopHat.

            #Calculate the stiffness of a single screw-fastened connection.
            Ka, ψ, α, β, Ke = Connections.cfs_trans_screwfastened_k(t_roof_deck, t_top_hat, E_roof_deck, E_top_hat, Fss, Fu_roof_deck, Fu_top_hat, deck_top_hat_fastener_diameter)

            #Convert the discrete stiffness to a distributed stiffness, divide by the fastener spacing.
            kx = Ke / deck_top_hat_fastener_spacing

            #Collect all the outputs.
            bracing_data[i] = PurlinLine.BracingData(kp, kϕ, kϕ_dist, kx, Lcrd, Lm)

        end

    elseif top_hat_purlin_line.inputs.new_deck_details[1] == "MR-24"  
        
        #There is no deck pullout stiffness needed here.
        kp = 0.0

        #Define the standing seam roof clip spacing.
        standing_seam_clip_spacing = top_hat_purlin_line.inputs.new_deck_details[2]

        #Define the standing seam roof clip height.
        standing_seam_clip_height = top_hat_purlin_line.inputs.new_deck_details[3]

        #Define the distance between clips as the distortional discrete bracing length.
        Lm = standing_seam_clip_spacing

        if standing_seam_clip_height == 2.25

            kϕ_standing_seam = 0.200  #From Seek et al. 2021, short floating clip, not exactly MR-24, kip-in/rad/in, https://www.researchgate.net/publication/349693825_Effective_standoff_in_standing_seam_roof_systems
            kx_standing_seam = 0.002  #From Cronin and Moen (2012), Figure 4.8  kips/in/in, https://vtechworks.lib.vt.edu/bitstream/handle/10919/18711/Flexural%20Capacity%20Prediction%20Method%20for%20an%20Open%20Web%20Joist%20Laterally%20Braced%20by%20a%20Standing%20Seam%20Roof%20System%20R10.pdf?sequence=1&isAllowed=y
        
        end

        #Loop over all the purlin segments in the line.  
        #Assume that there is a RoofHugger segment for every purlin segment.
        for i = 1:num_purlin_segments

            #Define the section property index associated with purlin segment i.
            section_index = top_hat_purlin_line.inputs.segments[i][2]

            #Define the material property index associated with purlin segment i.
            material_index = top_hat_purlin_line.inputs.segments[i][3]

            #Define the standing seam roof distributed clip stiffness.
            kϕ = kϕ_standing_seam

            #Define RoofHugger steel Poisson's ratio.
            μ_top_hat = top_hat_purlin_line.inputs.top_hat_material_properties[material_index][2]

            #Define the RoofHugger top flange width.
            b_top = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[section_index][3]

            #Define RoofHugger base metal thickness.
            t_top_hat = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[section_index][1]

            #Define out-to-out RoofHugger web depth.
            ho = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[section_index][4]

            #Define RoofHugger top flange lip length.
            d_top = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[section_index][2]

            #Define RoofHugger top flange lip angle from the horizon, in degrees.
            θ_top = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[section_index][9]

            #Apply Cee or Zee binary.   Assume the Roof Hugger behaves like a Z for this stiffness calculation.
            CorZ = 1

            #Calculate the RoofHugger distortional buckling half-wavelength.

            #Calculate top flange + lip section properties.
            Af, Jf, Ixf, Iyf, Ixyf, Cwf, xof,  hxf, hyf, yof = AISIS10016.table23131(CorZ, t_top_hat, b_top, d_top, θ_top)

            #Calculate the RoofHugger distortional buckling half-wavelength.
            Lcrd, L = AISIS10016.app23334(ho, μ_top_hat, t_top_hat, Ixf, xof, hxf, Cwf, Ixyf, Iyf, Lm)

            #If Lcrd is longer than the fastener spacing, then the distortional buckling will be restrained by the deck.
            if Lcrd >= Lm
                kϕ_dist = kϕ
            else
                kϕ_dist = 0.0
            end

            #Define standing seam roof lateral stiffness.
            kx = kx_standing_seam

            #Collect all the outputs.
            bracing_data[i] = PurlinLine.BracingData(kp, kϕ, kϕ_dist, kx, Lcrd, Lm)

        end

    elseif top_hat_purlin_line.inputs.new_deck_details[1] == "vertical leg standing seam"  
        
        #There is no deck pullout stiffness needed here.
        kp = 0.0

        #Define the standing seam roof clip spacing.
        standing_seam_clip_spacing = top_hat_purlin_line.inputs.new_deck_details[2]

        #Define the standing seam roof clip height.
        standing_seam_clip_height = top_hat_purlin_line.inputs.new_deck_details[3]

        #Define the distance between clips as the distortional discrete bracing length.
        Lm = standing_seam_clip_spacing

        if standing_seam_clip_height == 2.25

            kϕ_standing_seam = 0.200  #From Seek et al. 2021, short floating clip, not exactly MR-24, kip-in/rad/in, https://www.researchgate.net/publication/349693825_Effective_standoff_in_standing_seam_roof_systems
            kx_standing_seam = 0.002  #From Cronin and Moen (2012), Figure 4.8  kips/in/in, https://vtechworks.lib.vt.edu/bitstream/handle/10919/18711/Flexural%20Capacity%20Prediction%20Method%20for%20an%20Open%20Web%20Joist%20Laterally%20Braced%20by%20a%20Standing%20Seam%20Roof%20System%20R10.pdf?sequence=1&isAllowed=y
        
        end

        #Loop over all the purlin segments in the line.  
        #Assume that there is a RoofHugger segment for every purlin segment.
        for i = 1:num_purlin_segments

            #Define the section property index associated with purlin segment i.
            section_index = top_hat_purlin_line.inputs.segments[i][2]

            #Define the material property index associated with purlin segment i.
            material_index = top_hat_purlin_line.inputs.segments[i][3]

            #Define the standing seam roof distributed clip stiffness.
            kϕ = kϕ_standing_seam

            #Define RoofHugger steel Poisson's ratio.
            μ_top_hat = top_hat_purlin_line.inputs.top_hat_material_properties[material_index][2]

            #Define the RoofHugger top flange width.
            b_top = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[section_index][3]

            #Define RoofHugger base metal thickness.
            t_top_hat = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[section_index][1]

            #Define out-to-out RoofHugger web depth.
            ho = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[section_index][4]

            #Define RoofHugger top flange lip length.
            d_top = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[section_index][2]

            #Define RoofHugger top flange lip angle from the horizon, in degrees.
            θ_top = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[section_index][9]

            #Apply Cee or Zee binary.   Assume the Roof Hugger behaves like a Z for this stiffness calculation.
            CorZ = 1

            #Calculate the RoofHugger distortional buckling half-wavelength.

            #Calculate top flange + lip section properties.
            Af, Jf, Ixf, Iyf, Ixyf, Cwf, xof,  hxf, hyf, yof = AISIS10016.table23131(CorZ, t_top_hat, b_top, d_top, θ_top)

            #Calculate the RoofHugger distortional buckling half-wavelength.
            Lcrd, L = AISIS10016.app23334(ho, μ_top_hat, t_top_hat, Ixf, xof, hxf, Cwf, Ixyf, Iyf, Lm)

            #If Lcrd is longer than the fastener spacing, then the distortional buckling will be restrained by the deck.
            if Lcrd >= Lm
                kϕ_dist = kϕ
            else
                kϕ_dist = 0.0
            end

            #Define standing seam roof lateral stiffness.
            kx = kx_standing_seam

            #Collect all the outputs.
            bracing_data[i] = PurlinLine.BracingData(kp, kϕ, kϕ_dist, kx, Lcrd, Lm)

        end

    elseif top_hat_purlin_line.inputs.new_deck_details[1] == "no deck"
              
        #Loop over all the purlin segments in the line.  
        #Assume that there is a TopHat segment for every purlin segment.
        for i = 1:num_purlin_segments

            #Define the section property index associated with purlin segment i.
            section_index = top_hat_purlin_line.inputs.segments[i][2]

            #Define the material property index associated with purlin segment i.
            material_index = top_hat_purlin_line.inputs.segments[i][3]

            #Define TopHat steel elastic modulus.
            E_top_hat = top_hat_purlin_line.inputs.top_hat_material_properties[material_index][1]

            #Define TopHat steel Poisson's ratio.
            μ_top_hat = top_hat_purlin_line.inputs.top_hat_material_properties[material_index][2]

            #Define TopHat steel ultimate stress.
            Fu_top_hat = top_hat_purlin_line.inputs.top_hat_material_properties[material_index][4]

            #Define the TopHat top flange width.
            b_top = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[section_index][3]

            #Define TopHat base metal thickness.
            t_top_hat = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[section_index][1]

            #Define out-to-out TopHat web depth.
            ho = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[section_index][4]

            #Define TopHat top flange lip length.
            d_top = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[section_index][2]

            #Define TopHat top flange lip angle from the horizon, in degrees.
            θ_top = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[section_index][9]

            #Apply Cee or Zee binary.   Assume the TopHat behaves like a Z for this stiffness calculation.
            CorZ = 1

            #Calculate top flange + lip section properties.
            Af, Jf, Ixf, Iyf, Ixyf, Cwf, xof, hxf, hyf, yof = AISIS10016.table23131(CorZ, t_top_hat, b_top, d_top, θ_top)

            #Define the distance between fasteners as the distortional discrete bracing length.  There is no deck or fasteners in this case, so set Lm = length of purlin line.
            num_segments = size(top_hat_purlin_line.inputs.segments)[1]
            Lm = sum([top_hat_purlin_line.inputs.segments[i][1] for i = 1:num_segments])

            #Calculate the TopHat distortional buckling half-wavelength.
            Lcrd, L = AISIS10016.app23334(ho, μ_top_hat, t_top_hat, Ixf, xof, hxf, Cwf, Ixyf, Iyf, Lm)

            #Collect all the outputs.
            bracing_data[i] = PurlinLine.BracingData(0.0, 0.0, 0.0, 0.0, Lcrd, Lm)

        end

    end

    return bracing_data

end



function calculate_elastic_buckling_properties(top_hat_purlin_line)

    num_purlin_segments = size(top_hat_purlin_line.inputs.segments)[1]

    #Initialize vectors that will carry output.
    local_buckling_xx_pos = Array{PurlinLine.ElasticBucklingData, 1}(undef, num_purlin_segments)
    local_buckling_xx_neg = Array{PurlinLine.ElasticBucklingData, 1}(undef, num_purlin_segments)
    
    local_buckling_yy_pos = Array{PurlinLine.ElasticBucklingData, 1}(undef, num_purlin_segments)
    local_buckling_yy_neg = Array{PurlinLine.ElasticBucklingData, 1}(undef, num_purlin_segments)
    
    distortional_buckling_xx_pos = Array{PurlinLine.ElasticBucklingData, 1}(undef, num_purlin_segments)
    distortional_buckling_xx_neg = Array{PurlinLine.ElasticBucklingData, 1}(undef, num_purlin_segments)
 
    #Loop over all the purlin segments in the line.
    for i = 1:num_purlin_segments

        #Define the section property index associated with purlin segment i.
        section_index = top_hat_purlin_line.inputs.segments[i][2]

        #Define the material property index associated with purlin segment i.
        material_index = top_hat_purlin_line.inputs.segments[i][3]
        
        #Map section properties to CUFSM.
        A = top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].section_properties.A
        xcg = top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].section_properties.xc
        zcg = top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].section_properties.yc
        Ixx = top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].section_properties.Ixx
        Izz = top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].section_properties.Iyy
        Ixz = top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].section_properties.Ixy
        thetap = rad2deg(top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].section_properties.θ)
        I11 = top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].section_properties.I1
        I22 = top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].section_properties.I2
        unsymm = 0  #Sets Ixz=0 if unsymm = 0

        #Define the number of cross-section nodes.
        num_cross_section_nodes = size(top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].node_geometry)[1]

        #Initialize CUFSM node matrix.
        node = zeros(Float64, (num_cross_section_nodes, 8))

        #Add node numbers to node matrix.
        node[:, 1] .= 1:num_cross_section_nodes

        #Add nodal coordinates to node matrix.
        node[:, 2:3] .= top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].node_geometry

        #Add nodal restraints to node matrix.
        node[:, 4:7] .= ones(num_cross_section_nodes,4)

        #Define number of cross-section elements.
        num_cross_section_elements = size(top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].element_definitions)[1]

        #Initialize CUFSM elem matrix.
        elem = zeros(Float64, (num_cross_section_elements, 5))

        #Add element numbers to elem matrix.
        elem[:, 1] = 1:num_cross_section_elements

        #Add element connectivity and thickness to elem matrix.
        elem[:, 2:4] .= top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].element_definitions

        #Add element material reference to elem matrix.
        elem[:, 5] .= ones(num_cross_section_elements) * 100
                                
        #Find the purlin top flange centerline node.
        #lip curve bottom_flange curve web curve top_flange
        center_top_flange_purlin_node = sum(top_hat_purlin_line.purlin_cross_section_data[section_index].n[1:3]) + sum(top_hat_purlin_line.purlin_cross_section_data[section_index].n_radius[1:3]) + floor(Int, top_hat_purlin_line.purlin_cross_section_data[section_index].n[4]/2) + 1  #This floor command is a little dangerous.

        #Find the TopHat top flange centerline nodes.
        num_purlin_nodes = size(top_hat_purlin_line.purlin_cross_section_data[section_index].node_geometry)[1]
        center_minus_y_top_hat_flange_node = num_purlin_nodes + top_hat_purlin_line.top_hat_cross_section_data[section_index].n[1] + top_hat_purlin_line.top_hat_cross_section_data[section_index].n_radius[1] + floor(Int, top_hat_purlin_line.top_hat_cross_section_data[section_index].n[2] / 2) + 1
        center_plus_y_top_hat_flange_node = num_purlin_nodes + sum(top_hat_purlin_line.top_hat_cross_section_data[section_index].n[1:5]) + sum(top_hat_purlin_line.top_hat_cross_section_data[section_index].n_radius[1:5]) + floor(Int, top_hat_purlin_line.top_hat_cross_section_data[section_index].n[6] / 2) + 1

        #Set up springs in CUFSM.  There can be translational and rotational springs at the purlin top flange, and at each of the TopHat top flanges.
        springs = [1 center_top_flange_purlin_node 0 top_hat_purlin_line.bracing_data[i].kx 0 0 top_hat_purlin_line.bracing_data[i].kϕ_dist 0 0 0
                   2 center_minus_y_top_hat_flange_node 0 top_hat_purlin_line.new_deck_bracing_data[i].kx 0 0 top_hat_purlin_line.new_deck_bracing_data[i].kϕ_dist 0 0 0
                   3 center_plus_y_top_hat_flange_node 0 top_hat_purlin_line.new_deck_bracing_data[i].kx 0 0 top_hat_purlin_line.new_deck_bracing_data[i].kϕ_dist 0 0 0]
        
        #Constrain the TopHat bottom flange to the purlin top flange in all dof (x, z, y, and q).
        top_hat_bottom_flange_centerline_node = num_purlin_nodes + sum(top_hat_purlin_line.top_hat_cross_section_data[section_index].n[1:3]) + sum(top_hat_purlin_line.top_hat_cross_section_data[section_index].n_radius[1:3]) + floor(Int, top_hat_purlin_line.top_hat_cross_section_data[section_index].n[4] / 2) + 1
        
        #node#e DOFe coeff node#k DOFk
        constraints = [center_top_flange_purlin_node 1 1.0 top_hat_bottom_flange_centerline_node 1
                       center_top_flange_purlin_node 2 1.0 top_hat_bottom_flange_centerline_node 2
                       center_top_flange_purlin_node 3 1.0 top_hat_bottom_flange_centerline_node 3
                       center_top_flange_purlin_node 4 1.0 top_hat_bottom_flange_centerline_node 4]

        # constraints = 0

        #Assume here that purlin and TopHat have the same elastic modulus.
        E = top_hat_purlin_line.inputs.purlin_material_properties[material_index][1]
        ν = top_hat_purlin_line.inputs.purlin_material_properties[material_index][2]
        G = E / (2 *(1 + ν))
        prop = [100 E E ν ν G]

        neigs = 1  #just need the first mode 

        ###Local buckling - xx axis, positive 

        #Add reference stress to node matrix.

        #Define reference loads.  
        P = 0.0
        Mxx = 1.0  #assume centroidal moment always for now
        Mzz = 0.0
        M11 = 0.0
        M22 = 0.0

        #Define the TopHat web depth.
        h = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[section_index][4]  #this is a little dangerous
        length_inc = 5
        lengths = collect(0.25*h:0.75*h/length_inc:1.0*h)   #define to catch the local minimum

        CUFSM_local_xx_pos_data, Mcrℓ_xx_pos, Lcrℓ_xx_pos = PurlinLine.get_elastic_buckling(prop, deepcopy(node), elem, lengths, springs, constraints, neigs, P,Mxx,Mzz,M11,M22,A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,unsymm)   
        
        # mode_index = 4
        # shapes = CUFSM_local_xx_pos_data.shapes
        # scale_x = 1.0
        # scale_y = 1.0
        # # CUFSM.view_multi_branch_section_mode_shape(node, elem, shapes, mode_index, scale_x, scale_y)

        # half_wavelength = [curve[i,1][1] for i=1:length(lengths)]
        # load_factor = [curve[i,1][2] for i=1:length(lengths)]


        #Needed this deepcopy here to make struct work correctly.  Otherwise 'node' just kept changing.

        local_buckling_xx_pos[i] = PurlinLine.ElasticBucklingData(CUFSM_local_xx_pos_data, Lcrℓ_xx_pos, Mcrℓ_xx_pos)

        ###Local buckling - xx axis, negative 

        #Add reference stress to node matrix.

        #Define reference loads.  
        P = 0.0
        Mxx = -1.0  #assume centroidal moment always for now
        Mzz = 0.0
        M11 = 0.0
        M22 = 0.0

        #Use purlin web depth here.
        h = top_hat_purlin_line.inputs.purlin_cross_section_dimensions[section_index][5]  #this is a little dangerous
        length_inc = 5
        lengths = collect(0.25*h:0.75*h/length_inc:1.0*h)   #define to catch the local minimum

        CUFSM_local_xx_neg_data, Mcrℓ_xx_neg, Lcrℓ_xx_neg = PurlinLine.get_elastic_buckling(prop, deepcopy(node), elem, lengths, springs, constraints, neigs, P,Mxx,Mzz,M11,M22,A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,unsymm)

        local_buckling_xx_neg[i] = PurlinLine.ElasticBucklingData(CUFSM_local_xx_neg_data, Lcrℓ_xx_neg, Mcrℓ_xx_neg)


        ###local buckling - yy axis, positive  (compression at tip of purlin bottom flange lip)
        
        #Define reference loads.  
        P = 0.0
        Mxx = 0.0  
        Mzz = 1.0  #assume centroidal moment always for now
        M11 = 0.0
        M22 = 0.0

        #Try Lcrd of the purlin as a guide for finding the half-wavelength of the flange and lip (unstiffened element).
        #It turns out that for the case I studied, the purlin web buckled, not the purlin bottom flange lip.  Makes sense I guess.
        length_inc = 5
        lengths = collect(0.25 * top_hat_purlin_line.bracing_data[i].Lcrd:(1.0 * top_hat_purlin_line.bracing_data[i].Lcrd)/length_inc:1.25 * top_hat_purlin_line.bracing_data[i].Lcrd)   #define to catch the local minimum

        CUFSM_local_yy_pos_data, Mcrℓ_yy_pos, Lcrℓ_yy_pos = PurlinLine.get_elastic_buckling(prop, deepcopy(node), elem, lengths, springs, constraints, neigs, P,Mxx,Mzz,M11,M22,A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,unsymm)

    #                 mode_index = 3
    # shapes = CUFSM_local_yy_pos_data.shapes
    # scale_x = 1.0
    # scale_y = 1.0
    # CUFSM.view_multi_branch_section_mode_shape(node, elem, shapes, mode_index, scale_x, scale_y)


        local_buckling_yy_pos[i] = PurlinLine.ElasticBucklingData(CUFSM_local_yy_pos_data, Lcrℓ_yy_pos, Mcrℓ_yy_pos)
  
        ###local buckling - yy axis, negative  (compression in TopHat top flange lip)
        
        #Define reference loads.  
        P = 0.0
        Mxx = 0.0  
        Mzz = -1.0  #assume centroidal moment always for now
        M11 = 0.0
        M22 = 0.0
    
        length_inc = 5
        #Try Lcrd of the TopHat as a guide for finding the half-wavelength of the flange and lip (unstiffened element).
        lengths = collect(0.25 * top_hat_purlin_line.new_deck_bracing_data[i].Lcrd:(1.0 * top_hat_purlin_line.new_deck_bracing_data[i].Lcrd)/length_inc:1.25 * top_hat_purlin_line.new_deck_bracing_data[i].Lcrd)   #define to catch the local minimum

        CUFSM_local_yy_neg_data, Mcrℓ_yy_neg, Lcrℓ_yy_neg = PurlinLine.get_elastic_buckling(prop, deepcopy(node), elem, lengths, springs, constraints, neigs, P,Mxx,Mzz,M11,M22,A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,unsymm)

        local_buckling_yy_neg[i] = PurlinLine.ElasticBucklingData(CUFSM_local_yy_neg_data, Lcrℓ_yy_neg, Mcrℓ_yy_neg)

        ###Distortional buckling - xx axis, positive

        #Define reference loads.  
        P = 0.0
        Mxx = 1.0  #assume centroidal moment always for now
        Mzz = 0.0
        M11 = 0.0
        M22 = 0.0

        length_inc = 5
        #Use TopHat Lcrd here.
        lengths = collect(0.75 * top_hat_purlin_line.new_deck_bracing_data[i].Lcrd:(0.50 * top_hat_purlin_line.new_deck_bracing_data[i].Lcrd)/length_inc:1.25 * top_hat_purlin_line.new_deck_bracing_data[i].Lcrd)  #define to catch distortional minimum

        CUFSM_dist_pos_data, Mcrd_pos, Lcrd_pos_CUFSM = PurlinLine.get_elastic_buckling(prop, deepcopy(node), elem, lengths, springs, constraints, neigs, P,Mxx,Mzz,M11,M22,A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,unsymm)

        distortional_buckling_xx_pos[i] = PurlinLine.ElasticBucklingData(CUFSM_dist_pos_data, Lcrd_pos_CUFSM, Mcrd_pos)


 

         ###Distortional buckling - xx axis, negative

        #Define reference loads.  
        P = 0.0
        Mxx = -1.0  #assume centroidal moment always for now
        Mzz = 0.0
        M11 = 0.0
        M22 = 0.0

        length_inc = 5
        #Use purlin Lcrd here.
        lengths = collect(0.75 * top_hat_purlin_line.bracing_data[i].Lcrd:(0.50 * top_hat_purlin_line.bracing_data[i].Lcrd)/length_inc:1.25 * top_hat_purlin_line.bracing_data[i].Lcrd)  #define to catch distortional minimum

        CUFSM_dist_neg_data, Mcrd_neg, Lcrd_neg_CUFSM = PurlinLine.get_elastic_buckling(prop, deepcopy(node), elem, lengths, springs, constraints, neigs, P,Mxx,Mzz,M11,M22,A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,unsymm)

        distortional_buckling_xx_neg[i] = PurlinLine.ElasticBucklingData(CUFSM_dist_neg_data, Lcrd_neg_CUFSM, Mcrd_neg)

    end

    return local_buckling_xx_pos, local_buckling_xx_neg, local_buckling_yy_pos, local_buckling_yy_neg, distortional_buckling_xx_pos, distortional_buckling_xx_neg


end


function generate_top_hat_net_section_purlin_geometry(top_hat_purlin_cross_section_data, top_hat_cross_section_data, purlin_cross_section_data, purlin_cross_section_dimensions, top_hat_punch_out_dimensions)


    num_purlin_nodes = size(purlin_cross_section_data.node_geometry)[1]

    #Find all cross-section nodes that are in the TopHat punchout region.
    top_hat_node_geometry = top_hat_cross_section_data.node_geometry[:,1:2]
    hole_index_y = findall(x->x<=top_hat_punch_out_dimensions[2], top_hat_node_geometry[:,2])

    #Find all the nodes to the left and right of the TopHat webs.
    n_top_hat_web_minus = top_hat_cross_section_data.n[1] + top_hat_cross_section_data.n[2] + top_hat_cross_section_data.n_radius[1] + top_hat_cross_section_data.n_radius[2] + floor(Int, top_hat_cross_section_data.n[3]/2)

    n_top_hat_web_plus = sum(top_hat_cross_section_data.n[1:4]) + sum(top_hat_cross_section_data.n_radius[1:4]) + floor(Int, top_hat_cross_section_data.n[5]/2)

    top_hat_web_x_minus_location = top_hat_node_geometry[n_top_hat_web_minus, 1]
    top_hat_web_x_plus_location = top_hat_node_geometry[n_top_hat_web_plus, 1]
    hole_index_x = findall(x->(x<=top_hat_web_x_plus_location) & ((x>=top_hat_web_x_minus_location)), top_hat_node_geometry[:,1])

    #These are the nodes to be removed.
    hole_index = sort(intersect(hole_index_y, hole_index_x))

    #Shift the cross-section nodes to match up with the punch out dimensions.
    h_purlin = purlin_cross_section_dimensions[5]
    top_hat_purlin_node_geometry = top_hat_purlin_cross_section_data.node_geometry[:,1:2]
    top_hat_purlin_node_geometry[[hole_index[1] + num_purlin_nodes, hole_index[end] + num_purlin_nodes], 2] .=  h_purlin + top_hat_punch_out_dimensions[2]

    #Delete the nodes within the punchout region.
    remove_hole_nodes_index = hole_index[2:(end-1)] .+ num_purlin_nodes
    top_hat_purlin_node_geometry = top_hat_purlin_node_geometry[setdiff(1:end, remove_hole_nodes_index), :]

    #Remove elements at punchouts.
    top_hat_purlin_element_definitions = top_hat_purlin_cross_section_data.element_definitions
    remove_hole_elements_index = hole_index[1:end-1] .+ num_purlin_nodes .- 1
    top_hat_purlin_element_definitions = top_hat_purlin_element_definitions[setdiff(1:end, remove_hole_elements_index), :]

    #Update nodal connectivity.
    update_index = remove_hole_elements_index[1]
    num_removed_elements = length(remove_hole_elements_index)

    top_hat_purlin_element_definitions[update_index:end, 1:2] = top_hat_purlin_element_definitions[update_index:end, 1:2] .- num_removed_elements .+ 1
    
    return top_hat_purlin_node_geometry, top_hat_purlin_element_definitions

end


function define_top_hat_purlin_net_section(purlin_cross_section_dimensions, top_hat_cross_section_data, top_hat_plastic_cross_section_data, purlin_cross_section_data, purlin_plastic_cross_section_data, top_hat_purlin_cross_section_data, top_hat_purlin_plastic_cross_section_data, top_hat_punch_out_dimensions)

    num_purlin_sections = size(purlin_cross_section_dimensions)[1]

    cross_section_data = Vector{PurlinLine.CrossSectionData}(undef, num_purlin_sections)

    for i=1:num_purlin_sections

        top_hat_purlin_node_geometry, top_hat_purlin_element_definitions = generate_top_hat_net_section_purlin_geometry(top_hat_purlin_cross_section_data[i], top_hat_cross_section_data[i], purlin_cross_section_data[i], purlin_cross_section_dimensions[i],top_hat_punch_out_dimensions)

        #Calculate section properties of net section at a TopHat punchout.
        section_properties = CUFSM.cutwp_prop2(top_hat_purlin_node_geometry, top_hat_purlin_element_definitions)

        #Calculate TopHat+purlin plastic neutral axis and plastic modulus at the punchout net section.
        top_hat_purlin_plastic_node_geometry, top_hat_purlin_plastic_element_definitions = generate_top_hat_net_section_purlin_geometry(top_hat_purlin_plastic_cross_section_data[i], top_hat_plastic_cross_section_data[i], purlin_plastic_cross_section_data[i], purlin_cross_section_dimensions[i], top_hat_punch_out_dimensions)

        about_axis = "x"  #The strong axis plastic properties are needed for now.  
        top_hat_purlin_plastic_section_properties = StructuresKit.CrossSection.calculate_plastic_section_properties(top_hat_purlin_plastic_node_geometry, top_hat_purlin_plastic_element_definitions, about_axis)

        #Add cross section information to data structure.
        cross_section_data[i] = PurlinLine.CrossSectionData(top_hat_purlin_cross_section_data[i].n, top_hat_purlin_cross_section_data[i].n_radius, top_hat_purlin_node_geometry, top_hat_purlin_element_definitions, section_properties, top_hat_purlin_plastic_section_properties)

    end

    return cross_section_data

end



function calculate_net_section_local_buckling_properties(top_hat_purlin_line)
        
    num_purlin_segments = size(top_hat_purlin_line.inputs.segments)[1]

    #Initialize vectors that will carry output.
    local_buckling_xx_net_pos = Array{PurlinLine.ElasticBucklingData, 1}(undef, num_purlin_segments)
    
    #Loop over all the purlin segments in the line.
    for i = 1:num_purlin_segments

        #Define the section property index associated with purlin segment i.
        section_index = top_hat_purlin_line.inputs.segments[i][2]

        #Define the material property index associated with purlin segment i.
        material_index = top_hat_purlin_line.inputs.segments[i][3]
        
        #Map section properties to CUFSM.
        A = top_hat_purlin_line.top_hat_purlin_net_cross_section_data[section_index].section_properties.A
        xcg = top_hat_purlin_line.top_hat_purlin_net_cross_section_data[section_index].section_properties.xc
        zcg = top_hat_purlin_line.top_hat_purlin_net_cross_section_data[section_index].section_properties.yc
        Ixx = top_hat_purlin_line.top_hat_purlin_net_cross_section_data[section_index].section_properties.Ixx
        Izz = top_hat_purlin_line.top_hat_purlin_net_cross_section_data[section_index].section_properties.Iyy
        Ixz = top_hat_purlin_line.top_hat_purlin_net_cross_section_data[section_index].section_properties.Ixy
        thetap = rad2deg(top_hat_purlin_line.top_hat_purlin_net_cross_section_data[section_index].section_properties.θ)
        I11 = top_hat_purlin_line.top_hat_purlin_net_cross_section_data[section_index].section_properties.I1
        I22 = top_hat_purlin_line.top_hat_purlin_net_cross_section_data[section_index].section_properties.I2
        unsymm = 0  #Sets Ixz=0 if unsymm = 0

        #Define the number of cross-section nodes.
        num_cross_section_nodes = size(top_hat_purlin_line.top_hat_purlin_net_cross_section_data[section_index].node_geometry)[1]

        #Initialize CUFSM node matrix.
        node = zeros(Float64, (num_cross_section_nodes, 8))

        #Add node numbers to node matrix.
        node[:, 1] .= 1:num_cross_section_nodes

        #Add nodal coordinates to node matrix.
        node[:, 2:3] .= top_hat_purlin_line.top_hat_purlin_net_cross_section_data[section_index].node_geometry

        #Add nodal restraints to node matrix.
        node[:, 4:7] .= ones(num_cross_section_nodes,4)

        #Define number of cross-section elements.
        num_cross_section_elements = size(top_hat_purlin_line.top_hat_purlin_net_cross_section_data[section_index].element_definitions)[1]

        #Initialize CUFSM elem matrix.
        elem = zeros(Float64, (num_cross_section_elements, 5))

        #Add element numbers to elem matrix.
        elem[:, 1] = 1:num_cross_section_elements

        #Add element connectivity and thickness to elem matrix.
        elem[:, 2:4] .= top_hat_purlin_line.top_hat_purlin_net_cross_section_data[section_index].element_definitions

        #Add element material reference to elem matrix.
        elem[:, 5] .= ones(num_cross_section_elements) * 100
                                
        #Find the purlin top flange centerline node.
        #lip curve bottom_flange curve web curve top_flange
        center_top_flange_purlin_node =  sum(top_hat_purlin_line.purlin_cross_section_data[section_index].n[1:3]) + sum(top_hat_purlin_line.purlin_cross_section_data[section_index].n_radius[1:3]) + floor(Int, top_hat_purlin_line.purlin_cross_section_data[section_index].n[4]/2) + 1  #This floor command is a little dangerous.

        #Find the TopHat top flange centerline nodes.
        num_purlin_nodes = size(top_hat_purlin_line.purlin_cross_section_data[section_index].node_geometry)[1]
        center_minus_y_top_hat_flange_node = num_purlin_nodes + top_hat_purlin_line.top_hat_cross_section_data[section_index].n[1] + top_hat_purlin_line.top_hat_cross_section_data[section_index].n_radius[1] + floor(Int, top_hat_purlin_line.top_hat_cross_section_data[section_index].n[2] / 2) + 1
        
        num_top_hat_purlin_nodes = size(top_hat_purlin_line.top_hat_purlin_net_cross_section_data[section_index].node_geometry)[1]
        
        center_plus_y_top_hat_flange_node = num_top_hat_purlin_nodes -  top_hat_purlin_line.top_hat_purlin_net_cross_section_data[section_index].n[end] - top_hat_purlin_line.top_hat_purlin_net_cross_section_data[section_index].n_radius[end] - floor(Int, top_hat_purlin_line.top_hat_purlin_net_cross_section_data[section_index].n[end-1]/2)

        #Set up springs in CUFSM.  There can be translational and rotational springs at the purlin top flange, and at each of the TopHat top flanges.
        springs = [1 center_top_flange_purlin_node 0 top_hat_purlin_line.bracing_data[i].kx 0 0 top_hat_purlin_line.bracing_data[i].kϕ_dist 0 0 0
                2 center_minus_y_top_hat_flange_node 0 top_hat_purlin_line.new_deck_bracing_data[i].kx 0 0 top_hat_purlin_line.new_deck_bracing_data[i].kϕ_dist 0 0 0
                3 center_plus_y_top_hat_flange_node 0 top_hat_purlin_line.new_deck_bracing_data[i].kx 0 0 top_hat_purlin_line.new_deck_bracing_data[i].kϕ_dist 0 0 0]

        constraints = 0

        #Assume here that purlin and TopHat have the same elastic modulus.
        E = top_hat_purlin_line.inputs.purlin_material_properties[material_index][1]
        ν = top_hat_purlin_line.inputs.purlin_material_properties[material_index][2]
        G = E / (2 *(1 + ν))
        prop = [100 E E ν ν G]

        neigs = 1  #just need the first mode 

        ###Local buckling - xx axis, positive 

        #Add reference stress to node matrix.

        #Define reference loads.  
        P = 0.0
        Mxx = 1.0  #assume centroidal moment always for now
        Mzz = 0.0
        M11 = 0.0
        M22 = 0.0

        #Define the length of the punchout.
        L_hole = top_hat_purlin_line.inputs.top_hat_punch_out_dimensions[1]  

        #Assume the buckling half-wavelength is longer than TopHat punchout.  This seems reasonable based on wavelength studies where the mininum was around 15 in. for the net section TopHat + purlin model.  
        lengths = [L_hole]

        CUFSM_local_xx_net_pos_data, Mcrℓ_xx_net_pos, Lcrℓ_xx_net_pos = PurlinLine.get_elastic_buckling(prop, deepcopy(node), elem, lengths, springs, constraints, neigs, P,Mxx,Mzz,M11,M22,A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,unsymm)   

        #Needed this deepcopy here to make struct work correctly.  Otherwise 'node' just kept changing.

        local_buckling_xx_net_pos[i] = PurlinLine.ElasticBucklingData(CUFSM_local_xx_net_pos_data, Lcrℓ_xx_net_pos, Mcrℓ_xx_net_pos)

    end

    return local_buckling_xx_net_pos

end


function calculate_yielding_flexural_strength(top_hat_purlin_line)

    num_purlin_segments = size(top_hat_purlin_line.inputs.segments)[1]

    #Initialize a vectors that will hold all the outputs.
    yielding_flexural_strength_xx = Array{PurlinLine.YieldingFlexuralStrengthData, 1}(undef, num_purlin_segments)
    yielding_flexural_strength_xx_net = Array{PurlinLine.YieldingFlexuralStrengthData, 1}(undef, num_purlin_segments)
    yielding_flexural_strength_yy = Array{PurlinLine.YieldingFlexuralStrengthData, 1}(undef, num_purlin_segments)
    yielding_flexural_strength_free_flange_yy = Array{PurlinLine.YieldingFlexuralStrengthData, 1}(undef, num_purlin_segments)


    for i = 1:num_purlin_segments

        #Define the section property index associated with purlin segment i.
        section_index = top_hat_purlin_line.inputs.segments[i][2]

        #Define the material property index associated with purlin segment i.
        material_index = top_hat_purlin_line.inputs.segments[i][3]

        ###strong axis flexure, local-global interaction
        Fy_purlin = top_hat_purlin_line.inputs.purlin_material_properties[material_index][3]
        Fy_top_hat = top_hat_purlin_line.inputs.top_hat_material_properties[material_index][3]
        Ixx = top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].section_properties.Ixx
        cy_bottom = top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].section_properties.yc  #distance from neutral axis to bottom outer fiber

        top_hat_purlin_depth = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[section_index][4] + top_hat_purlin_line.inputs.purlin_cross_section_dimensions[section_index][5]
        cy_top = top_hat_purlin_depth - cy_bottom #distance from neutral axis to top outer fiber
        Sxx_pos = Ixx/cy_top
        Sxx_neg = Ixx/cy_bottom
        My_xx_pos = Fy_top_hat * Sxx_pos #TopHat at top fiber
        My_xx_neg = Fy_purlin * Sxx_neg  #purlin is at bottom fiber
        My_xx = minimum([My_xx_pos My_xx_neg])  #first yield criterion for AISI 

        yielding_flexural_strength_xx[i] = PurlinLine.YieldingFlexuralStrengthData(Sxx_pos, Sxx_neg, My_xx_pos, My_xx_neg, My_xx, 0.0)   #make eMy zero here since it is not used

        ###strong axis flexure, local-global interaction, net section at a punchout
        Ixx_net = top_hat_purlin_line.top_hat_purlin_net_cross_section_data[section_index].section_properties.Ixx
        cy_bottom_net = top_hat_purlin_line.top_hat_purlin_net_cross_section_data[section_index].section_properties.yc  #distance from neutral axis to bottom outer fiber

        cy_top_net = top_hat_purlin_depth - cy_bottom_net #distance from neutral axis to top outer fiber
        Sxx_pos_net = Ixx_net/cy_top_net
        Sxx_neg_net = Ixx_net/cy_bottom_net
        My_xx_pos_net = Fy_top_hat * Sxx_pos_net #TopHat at top fiber
        My_xx_neg_net = Fy_purlin * Sxx_neg_net  #purlin is at bottom fiber
        My_xx_net = minimum([My_xx_pos_net My_xx_neg_net])  #first yield criterion for AISI 

        yielding_flexural_strength_xx_net[i] = PurlinLine.YieldingFlexuralStrengthData(Sxx_pos_net, Sxx_neg_net, My_xx_pos_net, My_xx_neg_net, My_xx_net, 0.0)   #make eMy zero here since it is not used


        ###weak axis flexure, local-global interaction
        Iyy = top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].section_properties.Iyy

        #distance from neutral axis to (-x or left) outer fiber
        #Positive moment is applied when this outer fiber is compressed.
        cx_minusx = top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].section_properties.xc - minimum(top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].node_geometry[:,1])
        #distance from neutral axis to (+x or right) outer fiber
        #Negative moment is applied when this outer fiber is compressed.
        cx_plusx = maximum(top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].node_geometry[:,1]) - top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].section_properties.xc 
 
        Syy_pos = Iyy / cx_minusx
        Syy_neg = Iyy / cx_plusx
 
        My_yy_pos = Fy_purlin * Syy_pos  #outer fiber is the purlin
        My_yy_neg = Fy_top_hat *Syy_neg  #outer fiber is the TopHat.  This could be incorrect if TopHat yield stress is much higher than purlin yield stress!!!
        My_yy = minimum([My_yy_pos My_yy_neg])  #first yield criterion for AISI 
 
        yielding_flexural_strength_yy[i] = PurlinLine.YieldingFlexuralStrengthData(Syy_pos, Syy_neg, My_yy_pos, My_yy_neg, My_yy, 0.0)  #set eMy=0.0 for now

        ###free flange yy-axis, local-global interaction

        #define free flange properties
        Iyyf = top_hat_purlin_line.free_flange_cross_section_data[section_index].section_properties.Iyy

        #distance from neutral axis to (-x or left) outer fiber
        #Positive moment is applied when this outer fiber is compressed.
        cxf_minusx = top_hat_purlin_line.free_flange_cross_section_data[section_index].section_properties.xc - minimum(top_hat_purlin_line.free_flange_cross_section_data[section_index].node_geometry[:,1])
        #distance from neutral axis to (+x or right) outer fiber
        #Negative moment is applied when this outer fiber is compressed.
        cxf_plusx = maximum(top_hat_purlin_line.free_flange_cross_section_data[section_index].node_geometry[:,1]) - top_hat_purlin_line.free_flange_cross_section_data[section_index].section_properties.xc 

        Syy_pos_free_flange = Iyyf / cxf_minusx
        Syy_neg_free_flange = Iyyf / cxf_plusx

        My_yy_pos_free_flange = Fy_purlin * Syy_pos_free_flange
        My_yy_neg_free_flange = Fy_purlin * Syy_neg_free_flange
        My_yy_free_flange = minimum([My_yy_pos_free_flange My_yy_neg_free_flange])  #first yield criterion for AISI 

        #Factored yield moment is needed for the free flange to perform AISI interaction checks.

        if top_hat_purlin_line.inputs.design_code == "AISI S100-16 ASD"
            ASDorLRFD = 0
        elseif top_hat_purlin_line.inputs.design_code == "AISI S100-16 LRFD"
            ASDorLRFD = 1
        elseif top_hat_purlin_line.inputs.design_code == "AISI S100-16 nominal"
            ASDorLRFD = 2
        end

        Mcrℓ_yy_free_flange = 10.0^10 #Make this a big number so we just get back eMy
        My_yy_free_flange, eMy_yy_free_flange = AISIS10016.f321(My_yy_free_flange, Mcrℓ_yy_free_flange, ASDorLRFD)

        yielding_flexural_strength_free_flange_yy[i] = PurlinLine.YieldingFlexuralStrengthData(Syy_pos_free_flange, Syy_neg_free_flange, My_yy_pos_free_flange, My_yy_neg_free_flange, My_yy_free_flange, eMy_yy_free_flange)

    end

    return yielding_flexural_strength_xx, yielding_flexural_strength_xx_net, yielding_flexural_strength_yy, yielding_flexural_strength_free_flange_yy

end

function calculate_local_global_flexural_strength(top_hat_purlin_line)

    num_purlin_segments = size(top_hat_purlin_line.inputs.segments)[1]

    #Initialize a vectors that will hold all the outputs.
    local_global_flexural_strength_xx_no_hole = Array{PurlinLine.LocalGlobalFlexuralStrengthData, 1}(undef, num_purlin_segments)
    local_global_flexural_strength_xx_hole = Array{PurlinLine.LocalGlobalFlexuralStrengthData, 1}(undef, num_purlin_segments)
    local_global_flexural_strength_xx = Array{PurlinLine.LocalGlobalFlexuralStrengthData, 1}(undef, num_purlin_segments)
    local_global_flexural_strength_yy = Array{PurlinLine.LocalGlobalFlexuralStrengthData, 1}(undef, num_purlin_segments)
    local_global_flexural_strength_free_flange_yy = Array{PurlinLine.LocalGlobalFlexuralStrengthData, 1}(undef, num_purlin_segments)


    if top_hat_purlin_line.inputs.design_code == "AISI S100-16 ASD"
        ASDorLRFD = 0
    elseif top_hat_purlin_line.inputs.design_code == "AISI S100-16 LRFD"
        ASDorLRFD = 1
    elseif top_hat_purlin_line.inputs.design_code == "AISI S100-16 nominal"
        ASDorLRFD = 2
    end

    for i = 1:num_purlin_segments

        section_index = top_hat_purlin_line.inputs.segments[i][2]
        material_index = top_hat_purlin_line.inputs.segments[i][3]

        #Set Mne=My and handle global buckling with second order analysis.
        Mne_xx = top_hat_purlin_line.yielding_flexural_strength_xx[i].My  

        ###Work on positive flexure. Define Mnℓ as the mininum of [Mnℓ_no_hole, Mnℓ_hole] as suggested in AISI S100-16.
        
        #Start with section away from punchout.   Consider inelastic reserve.

        Mcrℓ_pos_no_hole = top_hat_purlin_line.local_buckling_xx_pos[i].Mcr

        λ_ℓ_pos_no_hole = sqrt(Mne_xx/Mcrℓ_pos_no_hole)

        if λ_ℓ_pos_no_hole < 0.776   #inelastic reserve is in play

            Sc = top_hat_purlin_line.yielding_flexural_strength_xx[i].S_pos
            St = top_hat_purlin_line.yielding_flexural_strength_xx[i].S_neg
            Z =  top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].plastic_section_properties.Z
            Fy_top_hat = top_hat_purlin_line.inputs.top_hat_material_properties[material_index][3]
            Fy_purlin = top_hat_purlin_line.inputs.purlin_material_properties[material_index][3]
            Fy = minimum([Fy_top_hat, Fy_purlin])  #take the minimum for now, can improve later
            lambda_l, Cyl, Mp, Myc, Myt3, Mnℓ_xx_pos_no_hole, eMnℓ_xx_pos_no_hole = AISIS10016.f323(Mne_xx, Mcrℓ_pos_no_hole, Sc, St, Z, Fy, ASDorLRFD)

        else

            Mnℓ_xx_pos_no_hole, eMnℓ_xx_pos_no_hole =  AISIS10016.f321(Mne_xx, Mcrℓ_pos_no_hole, ASDorLRFD)
        
        end

        #Now work on section at a punchout. Consider inelastic reserve.

        #Define Mnℓ as the mininum of [Mnℓ_no_hole, Mnℓ_hole] as suggested in AISI S100-16.
        My_net = top_hat_purlin_line.yielding_flexural_strength_xx_net[i].My
        Mcrℓ_pos_hole = top_hat_purlin_line.local_buckling_xx_net_pos[i].Mcr
    
        λ_ℓ_pos_hole = sqrt(Mne_xx/Mcrℓ_pos_hole)

        if λ_ℓ_pos_hole < 0.776   #inelastic reserve is in play

            Sc = top_hat_purlin_line.yielding_flexural_strength_xx_net[i].S_pos
            St = top_hat_purlin_line.yielding_flexural_strength_xx_net[i].S_neg
            Z =  top_hat_purlin_line.top_hat_purlin_net_cross_section_data[section_index].plastic_section_properties.Z
            Fy_top_hat = top_hat_purlin_line.inputs.top_hat_material_properties[material_index][3]
            Fy_purlin = top_hat_purlin_line.inputs.purlin_material_properties[material_index][3]
            Fy = minimum([Fy_top_hat, Fy_purlin])  #take the minimum for now, can improve later
            lambda_l, Cyl, Mp, Myc, Myt3, Mnℓ_xx_pos_hole, eMnℓ_xx_pos_hole = AISIS10016.f323(My_net, Mcrℓ_pos_hole, Sc, St, Z, Fy, ASDorLRFD)

        else

            Mnℓ_xx_pos_hole, eMnℓ_xx_pos_hole =  AISIS10016.f322(Mne_xx, Mcrℓ_pos_hole, My_net, ASDorLRFD)
        
        end

        Mnℓ_xx_pos = minimum([Mnℓ_xx_pos_no_hole, Mnℓ_xx_pos_hole])
        eMnℓ_xx_pos = minimum([eMnℓ_xx_pos_no_hole, eMnℓ_xx_pos_hole])

        ##### Now negative flexure.

        #Define Mnℓ as the mininum of [Mnℓ_no_hole, Mnℓ_hole] as suggested in AISI S100-16.
        Mcrℓ_neg_no_hole = top_hat_purlin_line.local_buckling_xx_neg[i].Mcr

        λ_ℓ_neg_no_hole = sqrt(Mne_xx/Mcrℓ_neg_no_hole)

        if λ_ℓ_neg_no_hole < 0.776   #inelastic reserve is in play

            Sc = top_hat_purlin_line.yielding_flexural_strength_xx[i].S_neg
            St = top_hat_purlin_line.yielding_flexural_strength_xx[i].S_pos
            Z =  top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].plastic_section_properties.Z
            Fy_top_hat = top_hat_purlin_line.inputs.top_hat_material_properties[material_index][3]
            Fy_purlin = top_hat_purlin_line.inputs.purlin_material_properties[material_index][3]
            Fy = minimum([Fy_top_hat, Fy_purlin])  #take the minimum for now, can improve later
            lambda_l, Cyl, Mp, Myc, Myt3, Mnℓ_xx_neg_no_hole, eMnℓ_xx_neg_no_hole = AISIS10016.f323(Mne_xx, Mcrℓ_neg_no_hole, Sc, St, Z, Fy, ASDorLRFD)

        else

            Mnℓ_xx_neg_no_hole, eMnℓ_xx_neg_no_hole =  AISIS10016.f321(Mne_xx, Mcrℓ_neg_no_hole, ASDorLRFD)
        
        end
        
        #At a TopHat punchout, net section yielding is the upper limit for negative flexure.  
        Mnℓ_xx_neg_hole, eMnℓ_xx_neg_hole =  AISIS10016.f322(Mne_xx, Mcrℓ_neg_no_hole, My_net, ASDorLRFD)

        #Find the minimum of the strengths at the punchout or away from the punchout.
        Mnℓ_xx_neg = minimum([Mnℓ_xx_neg_no_hole, Mnℓ_xx_neg_hole])
        eMnℓ_xx_neg = minimum([eMnℓ_xx_neg_no_hole, eMnℓ_xx_neg_hole])


        #Add no hole positive and negative flexural strength to data structure. 
        local_global_flexural_strength_xx_no_hole[i] = PurlinLine.LocalGlobalFlexuralStrengthData(Mne_xx, Mnℓ_xx_pos_no_hole, Mnℓ_xx_neg_no_hole, eMnℓ_xx_pos_no_hole, eMnℓ_xx_neg_no_hole)

        #Add net section positive and negative flexural strength to data structure. 
        local_global_flexural_strength_xx_hole[i] = PurlinLine.LocalGlobalFlexuralStrengthData(My_net, Mnℓ_xx_pos_hole, Mnℓ_xx_neg_hole, eMnℓ_xx_pos_hole, eMnℓ_xx_neg_hole)  

        #Add the governing (hole or no hole) positive and negative strong xx flexural strengths to the data structure.
        local_global_flexural_strength_xx[i] = PurlinLine.LocalGlobalFlexuralStrengthData(Mne_xx, Mnℓ_xx_pos, Mnℓ_xx_neg, eMnℓ_xx_pos, eMnℓ_xx_neg)

        ###weak axis flexure, local-global interaction
        Mne_yy = top_hat_purlin_line.yielding_flexural_strength_yy[i].My

        Mnℓ_yy_pos, eMnℓ_yy_pos = AISIS10016.f321(Mne_yy, top_hat_purlin_line.local_buckling_yy_pos[i].Mcr, ASDorLRFD)

        Mnℓ_yy_neg, eMnℓ_yy_neg = AISIS10016.f321(Mne_yy, top_hat_purlin_line.local_buckling_yy_neg[i].Mcr, ASDorLRFD)

        local_global_flexural_strength_yy[i] = PurlinLine.LocalGlobalFlexuralStrengthData(Mne_yy, Mnℓ_yy_pos, Mnℓ_yy_neg, eMnℓ_yy_pos, eMnℓ_yy_neg)


        ###free flange yy-axis, local-global interaction
        Mne_yy_free_flange = top_hat_purlin_line.yielding_flexural_strength_free_flange_yy[i].My 

        #Assume no local buckling for now in the free flange strength calculation.  Set Mcrℓ to Mne times a big number. 

        Mnℓ_yy_pos_free_flange, eMnℓ_yy_pos_free_flange = AISIS10016.f321(Mne_yy_free_flange, Mne_yy_free_flange * 1000, ASDorLRFD)

        Mnℓ_yy_neg_free_flange, eMnℓ_yy_neg_free_flange = AISIS10016.f321(Mne_yy_free_flange, Mne_yy_free_flange * 1000, ASDorLRFD)

        local_global_flexural_strength_free_flange_yy[i] = PurlinLine.LocalGlobalFlexuralStrengthData(Mne_yy_free_flange, Mnℓ_yy_pos_free_flange, Mnℓ_yy_neg_free_flange, eMnℓ_yy_pos_free_flange, eMnℓ_yy_neg_free_flange)

    end

    return local_global_flexural_strength_xx_no_hole, local_global_flexural_strength_xx_hole,local_global_flexural_strength_xx, local_global_flexural_strength_yy, local_global_flexural_strength_free_flange_yy

end



function define_top_hat_purlin_distortional_net_section(top_hat_purlin_line)

    num_purlin_segments = size(top_hat_purlin_line.inputs.segments)[1]

    cross_section_data = Vector{PurlinLine.CrossSectionData}(undef, num_purlin_segments)

    #Define length of TopHat punchout.
    L_hole = top_hat_purlin_line.inputs.top_hat_punch_out_dimensions[1]

    for i=1:num_purlin_segments

        #Define the section property index associated with purlin segment i.
        section_index = top_hat_purlin_line.inputs.segments[i][2]

        ##Find all cross-section nodes that are in the TopHat webs.
        top_hat_purlin_node_geometry = top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].node_geometry[:,1:2]

        #Define the x range for the TopHat webs.

        #Work from the purlin top flange centerline.
        center_top_flange_purlin_node =  sum(top_hat_purlin_line.purlin_cross_section_data[section_index].n[1:3]) + sum(top_hat_purlin_line.purlin_cross_section_data[section_index].n_radius[1:3]) + floor(Int, top_hat_purlin_line.purlin_cross_section_data[section_index].n[4]/2) + 1  

        #Define the TopHat bottom flange width.
        top_hat_bottom_flange_width = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[section_index][5]

        #Define the TopHat base metal thickness.
        t_top_hat = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[section_index][1]

        #Define the TopHat minus x web coordinate.
        top_hat_web_x_minus = top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].node_geometry[center_top_flange_purlin_node, 1] - top_hat_bottom_flange_width/2 + t_top_hat/2

        #Define the TopHat plus x web coordinate.
        top_hat_web_x_plus = top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].node_geometry[center_top_flange_purlin_node, 1] + top_hat_bottom_flange_width/2 - t_top_hat/2

        web_index_minus = findall(x->x≈top_hat_web_x_minus, top_hat_purlin_node_geometry[:,1])

        web_index_plus = findall(x->x≈top_hat_web_x_plus, top_hat_purlin_node_geometry[:,1])

        #Approximate the Lcrd for each cross-section.   
        Lcrd = top_hat_purlin_line.distortional_buckling_xx_pos[section_index].Lcr
        tr = AISIS10016.app2C2262(t_top_hat, L_hole, Lcrd)

        #Define element definitions.
        top_hat_purlin_element_definitions = top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].element_definitions

        #Define the element ranges where t becomes tr.
        web_element_index_minus = [web_index_minus[1] - 2; web_index_minus[1] - 1; web_index_minus[1:(end-1)]]

        web_element_index_plus = [web_index_plus[1] - 2; web_index_plus[1] - 1; web_index_plus[1:(end-1)]]

        #Update the element thicknesses in the TopHat.
        top_hat_purlin_element_definitions[web_element_index_minus, 3] .= tr
        top_hat_purlin_element_definitions[web_element_index_plus, 3] .= tr
    
        #Calculate section properties of section with reduced web thickness.
        section_properties = CUFSM.cutwp_prop2(top_hat_purlin_node_geometry, top_hat_purlin_element_definitions)
    
        #Add cross section information to data structure.
        cross_section_data[i] = PurlinLine.CrossSectionData(top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].n, top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].n_radius, top_hat_purlin_node_geometry, top_hat_purlin_element_definitions, section_properties, nothing)

    end

    return cross_section_data

end




function calculate_net_section_distortional_buckling_properties(top_hat_purlin_line)
        
    num_purlin_segments = size(top_hat_purlin_line.inputs.segments)[1]

    #Initialize vectors that will carry output.
    distortional_buckling_xx_net_pos = Array{PurlinLine.ElasticBucklingData, 1}(undef, num_purlin_segments)
    
    #Loop over all the purlin segments in the line.
    for i = 1:num_purlin_segments

        #Define the section property index associated with purlin segment i.
        section_index = top_hat_purlin_line.inputs.segments[i][2]

        #Define the material property index associated with purlin segment i.
        material_index = top_hat_purlin_line.inputs.segments[i][3]
        
        #Map section properties to CUFSM.
        A = top_hat_purlin_line.top_hat_purlin_distortional_net_cross_section_data[section_index].section_properties.A
        xcg = top_hat_purlin_line.top_hat_purlin_distortional_net_cross_section_data[section_index].section_properties.xc
        zcg = top_hat_purlin_line.top_hat_purlin_distortional_net_cross_section_data[section_index].section_properties.yc
        Ixx = top_hat_purlin_line.top_hat_purlin_distortional_net_cross_section_data[section_index].section_properties.Ixx
        Izz = top_hat_purlin_line.top_hat_purlin_distortional_net_cross_section_data[section_index].section_properties.Iyy
        Ixz = top_hat_purlin_line.top_hat_purlin_distortional_net_cross_section_data[section_index].section_properties.Ixy
        thetap = rad2deg(top_hat_purlin_line.top_hat_purlin_distortional_net_cross_section_data[section_index].section_properties.θ)
        I11 = top_hat_purlin_line.top_hat_purlin_distortional_net_cross_section_data[section_index].section_properties.I1
        I22 = top_hat_purlin_line.top_hat_purlin_distortional_net_cross_section_data[section_index].section_properties.I2
        unsymm = 0  #Sets Ixz=0 if unsymm = 0

        #Define the number of cross-section nodes.
        num_cross_section_nodes = size(top_hat_purlin_line.top_hat_purlin_distortional_net_cross_section_data[section_index].node_geometry)[1]

        #Initialize CUFSM node matrix.
        node = zeros(Float64, (num_cross_section_nodes, 8))

        #Add node numbers to node matrix.
        node[:, 1] .= 1:num_cross_section_nodes

        #Add nodal coordinates to node matrix.
        node[:, 2:3] .= top_hat_purlin_line.top_hat_purlin_distortional_net_cross_section_data[section_index].node_geometry

        #Add nodal restraints to node matrix.
        node[:, 4:7] .= ones(num_cross_section_nodes,4)

        #Define number of cross-section elements.
        num_cross_section_elements = size(top_hat_purlin_line.top_hat_purlin_distortional_net_cross_section_data[section_index].element_definitions)[1]

        #Initialize CUFSM elem matrix.
        elem = zeros(Float64, (num_cross_section_elements, 5))

        #Add element numbers to elem matrix.
        elem[:, 1] = 1:num_cross_section_elements

        #Add element connectivity and thickness to elem matrix.
        elem[:, 2:4] .= top_hat_purlin_line.top_hat_purlin_distortional_net_cross_section_data[section_index].element_definitions

        #Add element material reference to elem matrix.
        elem[:, 5] .= ones(num_cross_section_elements) * 100
                                
        #Find the purlin top flange centerline node.
        #lip curve bottom_flange curve web curve top_flange
        center_top_flange_purlin_node =  sum(top_hat_purlin_line.purlin_cross_section_data[section_index].n[1:3]) + sum(top_hat_purlin_line.purlin_cross_section_data[section_index].n_radius[1:3]) + floor(Int, top_hat_purlin_line.purlin_cross_section_data[section_index].n[4]/2) + 1  #This floor command is a little dangerous.
        
        #Find the TopHat top flange centerline nodes.
        num_purlin_nodes = size(top_hat_purlin_line.purlin_cross_section_data[section_index].node_geometry)[1]
        center_minus_y_top_hat_flange_node = num_purlin_nodes + top_hat_purlin_line.top_hat_cross_section_data[section_index].n[1] + top_hat_purlin_line.top_hat_cross_section_data[section_index].n_radius[1] + floor(Int, top_hat_purlin_line.top_hat_cross_section_data[section_index].n[2] / 2) + 1
        center_plus_y_top_hat_flange_node = num_purlin_nodes + sum(top_hat_purlin_line.top_hat_cross_section_data[section_index].n[1:5]) + sum(top_hat_purlin_line.top_hat_cross_section_data[section_index].n_radius[1:5]) + floor(Int, top_hat_purlin_line.top_hat_cross_section_data[section_index].n[6] / 2) + 1

        #Set up springs in CUFSM.  There can be translational and rotational springs at the purlin top flange, and at each of the TopHat top flanges.
        springs = [1 center_top_flange_purlin_node 0 top_hat_purlin_line.bracing_data[i].kx 0 0 top_hat_purlin_line.bracing_data[i].kϕ_dist 0 0 0
                   2 center_minus_y_top_hat_flange_node 0 top_hat_purlin_line.new_deck_bracing_data[i].kx 0 0 top_hat_purlin_line.new_deck_bracing_data[i].kϕ_dist 0 0 0
                   3 center_plus_y_top_hat_flange_node 0 top_hat_purlin_line.new_deck_bracing_data[i].kx 0 0 top_hat_purlin_line.new_deck_bracing_data[i].kϕ_dist 0 0 0]
        
        #Constrain the TopHat bottom flange to the purlin top flange in all dof (x, z, y, and q).
        top_hat_bottom_flange_centerline_node = num_purlin_nodes + sum(top_hat_purlin_line.top_hat_cross_section_data[section_index].n[1:3]) + sum(top_hat_purlin_line.top_hat_cross_section_data[section_index].n_radius[1:3]) + floor(Int, top_hat_purlin_line.top_hat_cross_section_data[section_index].n[4] / 2) + 1
        
        #node#e DOFe coeff node#k DOFk
        constraints = [center_top_flange_purlin_node 1 1.0 top_hat_bottom_flange_centerline_node 1
                       center_top_flange_purlin_node 2 1.0 top_hat_bottom_flange_centerline_node 2
                       center_top_flange_purlin_node 3 1.0 top_hat_bottom_flange_centerline_node 3
                       center_top_flange_purlin_node 4 1.0 top_hat_bottom_flange_centerline_node 4]

        #Assume here that purlin and TopHat have the same elastic modulus.
        E = top_hat_purlin_line.inputs.purlin_material_properties[material_index][1]
        ν = top_hat_purlin_line.inputs.purlin_material_properties[material_index][2]
        G = E / (2 *(1 + ν))
        prop = [100 E E ν ν G]

        neigs = 1  #just need the first mode 

        ###Local buckling - xx axis, positive 

        #Add reference stress to node matrix.

        #Define reference loads.  
        P = 0.0
        Mxx = 1.0  #assume centroidal moment always for now
        Mzz = 0.0
        M11 = 0.0
        M22 = 0.0

        #Define the distortional buckling half-wavelength.
        Lcrd = top_hat_purlin_line.distortional_buckling_xx_pos[i].Lcr 

        #Calculate the buckling load just at Lcrd.
        lengths = [Lcrd]

        CUFSM_distortional_xx_net_pos_data, Mcrd_xx_net_pos, Lcrd_xx_net_pos = PurlinLine.get_elastic_buckling(prop, deepcopy(node), elem, lengths, springs, constraints, neigs, P,Mxx,Mzz,M11,M22,A,xcg,zcg,Ixx,Izz,Ixz,thetap,I11,I22,unsymm)   

        #Needed this deepcopy here to make struct work correctly.  Otherwise 'node' just kept changing.

        distortional_buckling_xx_net_pos[i] = PurlinLine.ElasticBucklingData(CUFSM_distortional_xx_net_pos_data, Lcrd_xx_net_pos, Mcrd_xx_net_pos)

    end

    return distortional_buckling_xx_net_pos

end


function calculate_distortional_flexural_strength(top_hat_purlin_line)

    num_purlin_segments = size(top_hat_purlin_line.inputs.segments)[1]

    #Initialize a vectors that will hold all the outputs.
    distortional_flexural_strength_xx = Array{PurlinLine.DistortionalFlexuralStrengthData, 1}(undef, num_purlin_segments)

    if top_hat_purlin_line.inputs.design_code == "AISI S100-16 ASD"
        ASDorLRFD = 0
    elseif top_hat_purlin_line.inputs.design_code == "AISI S100-16 LRFD"
        ASDorLRFD = 1
    elseif top_hat_purlin_line.inputs.design_code == "AISI S100-16 nominal"
        ASDorLRFD = 2
    end

    for i = 1:num_purlin_segments

        Mnd_xx_pos, eMnd_xx_pos = AISIS10016.f411(top_hat_purlin_line.yielding_flexural_strength_xx[i].My, top_hat_purlin_line.distortional_buckling_xx_net_pos[i].Mcr, ASDorLRFD)  #use Mcrd_hole always here

        Mnd_xx_neg, eMnd_xx_neg = AISIS10016.f411(top_hat_purlin_line.yielding_flexural_strength_xx[i].My, top_hat_purlin_line.distortional_buckling_xx_neg[i].Mcr, ASDorLRFD)  #assume holes do not affect negative bending distortional buckling for TopHat + purlin

        distortional_flexural_strength_xx[i] = PurlinLine.DistortionalFlexuralStrengthData(Mnd_xx_pos, Mnd_xx_neg, eMnd_xx_pos, eMnd_xx_neg)

    end

    return distortional_flexural_strength_xx

end


function calculate_torsion_strength(top_hat_purlin_line)

    num_purlin_segments = size(top_hat_purlin_line.inputs.segments)[1]

    #Initialize a vector that will hold all the outputs.
    torsion_strength = Array{PurlinLine.TorsionStrengthData, 1}(undef, num_purlin_segments)
    
    if top_hat_purlin_line.inputs.design_code == "AISI S100-16 ASD"
        ASDorLRFD = 0
    elseif top_hat_purlin_line.inputs.design_code == "AISI S100-16 LRFD"
        ASDorLRFD = 1
    elseif top_hat_purlin_line.inputs.design_code == "AISI S100-16 nominal"
        ASDorLRFD = 2
    end

    for i = 1:num_purlin_segments

        #Define the section property index associated with purlin segment i.
        section_index = top_hat_purlin_line.inputs.segments[i][2]

        #Define the material property index associated with purlin segment i.
        material_index = top_hat_purlin_line.inputs.segments[i][3]
        
        Cw = top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].section_properties.Cw

        #This is the maximum magnitude of the warping stress function.  
        Wn = maximum(abs.(top_hat_purlin_line.top_hat_purlin_cross_section_data[section_index].section_properties.wn))

        Fy_purlin = top_hat_purlin_line.inputs.purlin_material_properties[material_index][3]
        Fy_top_hat = top_hat_purlin_line.inputs.top_hat_material_properties[material_index][3]
        Fy = minimum([Fy_purlin, Fy_top_hat])  #Maximum warping stress will be in the top right TopHat flange or the bottom purlin flange lip, so use the minimum yield stress here.

        Bn, eBn = AISIS10024.h411(Cw, Fy, Wn, ASDorLRFD)

        torsion_strength[i] = PurlinLine.TorsionStrengthData(Wn, Bn, eBn)

    end

    return torsion_strength

end


function calculate_shear_strength(top_hat_purlin_line)

    num_purlin_segments = size(top_hat_purlin_line.inputs.segments)[1]

    #Initialize a vector that will hold all the outputs.
    shear_strength_purlin = Array{PurlinLine.ShearStrengthData, 1}(undef, num_purlin_segments)
    shear_strength_top_hat = Array{PurlinLine.ShearStrengthData, 1}(undef, num_purlin_segments)
    shear_strength = Array{PurlinLine.ShearStrengthData, 1}(undef, num_purlin_segments)

    if top_hat_purlin_line.inputs.design_code == "AISI S100-16 ASD"
        ASDorLRFD = 0
    elseif top_hat_purlin_line.inputs.design_code == "AISI S100-16 LRFD"
        ASDorLRFD = 1
    elseif top_hat_purlin_line.inputs.design_code == "AISI S100-16 nominal"
        ASDorLRFD = 2
    end

    for i = 1:num_purlin_segments
    
        #Define the section property index associated with purlin segment i.
        section_index = top_hat_purlin_line.inputs.segments[i][2]

        #Define the material property index associated with purlin segment i.
        material_index = top_hat_purlin_line.inputs.segments[i][3]

        #Set a, the shear stiffener spacing, to the sum of the purlin segment lengths.  This assumes that shear stiffeners are not provided.
        sum_purlin_segments = sum([top_hat_purlin_line.inputs.segments[i][1] for i=1:size(top_hat_purlin_line.inputs.segments)[1]])
        a = sum_purlin_segments

        #Assume shear strength is Vn,purlin + Vn,top_hat.

        #Vn, purlin

        #Define base metal thickness.
        t_purlin = top_hat_purlin_line.inputs.purlin_cross_section_dimensions[section_index][2]

        #Define material properties.
        E_purlin = top_hat_purlin_line.inputs.purlin_material_properties[material_index][1]
        μ_purlin = top_hat_purlin_line.inputs.purlin_material_properties[material_index][2]
        Fy_purlin = top_hat_purlin_line.inputs.purlin_material_properties[material_index][3]

        #Depth of flat portion of web.
        full_web_depth_purlin = top_hat_purlin_line.inputs.purlin_cross_section_dimensions[section_index][5]
        bottom_flange_web_outside_radius_purlin = top_hat_purlin_line.inputs.purlin_cross_section_dimensions[section_index][14]
        top_flange_web_outside_radius_purlin = top_hat_purlin_line.inputs.purlin_cross_section_dimensions[section_index][15]
        h_flat_purlin = full_web_depth_purlin - bottom_flange_web_outside_radius_purlin - top_flange_web_outside_radius_purlin

        #Calculate plate buckling coefficient.
        kv_purlin  = AISIS10016.g233(a, h_flat_purlin)

        #Calculate shear buckling stress.
        Fcrv_purlin = AISIS10016.g232(E_purlin, μ_purlin, kv_purlin, h_flat_purlin, t_purlin)
        Vcr_purlin = AISIS10016.g231(h_flat_purlin, t_purlin, Fcrv_purlin)

        #Calculate shear buckling strength.
        Vn_purlin, eVn_purlin = AISIS10016.g21(h_flat_purlin, t_purlin, Fy_purlin, Vcr_purlin, ASDorLRFD)

        #Vn, TopHat

        #Define base metal thickness.
        t_top_hat = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[section_index][1]

        #Define material properties.
        E_top_hat = top_hat_purlin_line.inputs.top_hat_material_properties[material_index][1]
        μ_top_hat = top_hat_purlin_line.inputs.top_hat_material_properties[material_index][2]
        Fy_top_hat = top_hat_purlin_line.inputs.top_hat_material_properties[material_index][3]

        #Depth of flat portion of web.
        full_web_depth_top_hat = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[section_index][4]
        bottom_flange_web_outside_radius_top_hat = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[section_index][18]
        top_flange_web_outside_radius_top_hat = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[section_index][17]
        h_flat_top_hat = full_web_depth_top_hat - bottom_flange_web_outside_radius_top_hat - top_flange_web_outside_radius_top_hat

        #Calculate plate buckling coefficient.
        kv_top_hat  = AISIS10016.g233(a, h_flat_top_hat)

        #Calculate shear buckling stress.
        Fcrv_top_hat = AISIS10016.g232(E_top_hat, μ_top_hat, kv_top_hat, h_flat_top_hat, t_top_hat)
        Vcr_top_hat = AISIS10016.g231(h_flat_top_hat, t_top_hat, Fcrv_top_hat)

        #Calculate shear buckling strength for one web of TopHat.
        Vn_top_hat, eVn_top_hat = AISIS10016.g21(h_flat_top_hat, t_top_hat, Fy_top_hat, Vcr_top_hat, ASDorLRFD)

        Vn = Vn_purlin + 2 * Vn_top_hat
        eVn = eVn_purlin + 2 * eVn_top_hat

        shear_strength_purlin[i] = PurlinLine.ShearStrengthData(h_flat_purlin, kv_purlin, Fcrv_purlin, Vcr_purlin, Vn_purlin, eVn_purlin)
        shear_strength_top_hat[i] = PurlinLine.ShearStrengthData(h_flat_top_hat, kv_top_hat, Fcrv_top_hat, Vcr_top_hat, Vn_top_hat, eVn_top_hat)
        shear_strength[i] = PurlinLine.ShearStrengthData(0.0, 0.0, 0.0, 0.0, Vn, eVn)
    

    end

    return shear_strength_purlin, shear_strength_top_hat, shear_strength

end


function define(design_code, segments, spacing, roof_slope, purlin_cross_section_dimensions, top_hat_cross_section_dimensions, top_hat_punch_out_dimensions, purlin_material_properties, top_hat_material_properties, deck_details, deck_material_properties, new_deck_details, new_deck_material_properties, frame_flange_width, support_locations, purlin_frame_connections, bridging_locations)

    #Create the TopHatDesigner data structure.
    top_hat_purlin_line = TopHatDesignerObject()

    #Add TopHatDesigner user inputs to data structure.
    top_hat_purlin_line.inputs = TopHatDesigner.Inputs(design_code, segments, spacing, roof_slope, purlin_cross_section_dimensions, top_hat_cross_section_dimensions, top_hat_punch_out_dimensions, purlin_material_properties, top_hat_material_properties, new_deck_details, new_deck_material_properties, frame_flange_width, support_locations, purlin_frame_connections, bridging_locations)

    #Define TopHat cross-section data including nodal geometry, cross-section discretization and section properties.
    n = [4, 4, 6, 4, 6, 4, 4]
    n_radius = [4, 4, 4, 4, 4, 4]
    top_hat_purlin_line.top_hat_cross_section_data = define_top_hat_cross_sections(top_hat_purlin_line.inputs.top_hat_cross_section_dimensions, n, n_radius)

    #Create the PurlinLine data structure.
    purlin_line = PurlinLine.PurlinLineObject()

    #Capture PurlinLine inputs.
    purlin_line.inputs = PurlinLine.Inputs(design_code, segments, spacing, roof_slope, purlin_cross_section_dimensions, purlin_material_properties, deck_details, deck_material_properties, frame_flange_width, support_locations, purlin_frame_connections, bridging_locations)

    #Define the purlin cross-section discretization and calculate section properties.
    n = [4, 4, 5, 4, 4]
    n_radius = [4, 4, 4, 4]
    top_hat_purlin_line.purlin_cross_section_data = PurlinLine.define_purlin_section(purlin_line.inputs.cross_section_dimensions, n, n_radius)
    purlin_line.cross_section_data = top_hat_purlin_line.purlin_cross_section_data

    #Define the purlin free flange cross-section discretization and calculate section properties.
    n = [4, 4, 4]
    n_radius = [4, 4]
    top_hat_purlin_line.free_flange_cross_section_data = PurlinLine.define_purlin_free_flange_section(purlin_line.inputs.cross_section_dimensions, n, n_radius)
    purlin_line.free_flange_cross_section_data = top_hat_purlin_line.free_flange_cross_section_data

    #Define the purlin cross-section discretization for calculating plastic properties.
    n = [4, 4, 5, 4, 4] .* 10
    n_radius = [4, 4, 4, 4]
    purlin_plastic_cross_section_data = PurlinLine.define_purlin_section(purlin_line.inputs.cross_section_dimensions, n, n_radius)
   
    #Define TopHat cross-section discretizaton for calculating plastic properties.
    n = [4, 4, 6, 4, 6, 4, 4] .* 10
    n_radius = [4, 4, 4, 4, 4, 4]
    top_hat_plastic_cross_section_data = define_top_hat_cross_sections(top_hat_purlin_line.inputs.top_hat_cross_section_dimensions, n, n_radius)

    #Define TopHat + purlin properties.
    top_hat_purlin_line.top_hat_purlin_cross_section_data = define_top_hat_purlin_cross_sections(top_hat_purlin_line.inputs.purlin_cross_section_dimensions, top_hat_purlin_line.purlin_cross_section_data, top_hat_purlin_line.top_hat_cross_section_data, purlin_plastic_cross_section_data, top_hat_plastic_cross_section_data)
   
    #Define TopHat + purlin plastic discretization.
    top_hat_purlin_plastic_cross_section_data = define_top_hat_purlin_cross_sections(top_hat_purlin_line.inputs.purlin_cross_section_dimensions, purlin_plastic_cross_section_data, top_hat_plastic_cross_section_data, purlin_plastic_cross_section_data, top_hat_plastic_cross_section_data)

    #Calculate deck bracing properties.  This is for the purlin to deck.   
    #Assume for now that adding big screws and the TopHat does not change the stiffness properties. 
    top_hat_purlin_line.bracing_data = PurlinLine.define_deck_bracing_properties(purlin_line)
    purlin_line.bracing_data = PurlinLine.define_deck_bracing_properties(purlin_line)

    #Calculate free flange shear flow properties, including bracing stiffness from web and conversion factor from purlin line load to shear flow.
    #Assume shear flow is unchanged with addition of the TopHat.
    top_hat_purlin_line.free_flange_data = PurlinLine.calculate_free_flange_shear_flow_properties(purlin_line)


    #Calculate bracing properties for new roof deck to TopHat.
    top_hat_purlin_line.new_deck_bracing_data = define_new_deck_bracing_properties(top_hat_purlin_line)

    #Calculate elastic buckling properties for TopHat and purlin as a cross-section together.
    top_hat_purlin_line.local_buckling_xx_pos, top_hat_purlin_line.local_buckling_xx_neg, top_hat_purlin_line.local_buckling_yy_pos, top_hat_purlin_line.local_buckling_yy_neg, top_hat_purlin_line.distortional_buckling_xx_pos, top_hat_purlin_line.distortional_buckling_xx_neg = TopHatDesigner.calculate_elastic_buckling_properties(top_hat_purlin_line)

    #Calculate net section properties and add to TopHatDesigner data structure.
    top_hat_purlin_line.top_hat_purlin_net_cross_section_data = define_top_hat_purlin_net_section(top_hat_purlin_line.inputs.purlin_cross_section_dimensions, top_hat_purlin_line.top_hat_cross_section_data, top_hat_plastic_cross_section_data, top_hat_purlin_line.purlin_cross_section_data, purlin_plastic_cross_section_data,  top_hat_purlin_line.top_hat_purlin_cross_section_data, top_hat_purlin_plastic_cross_section_data, top_hat_purlin_line.inputs.top_hat_punch_out_dimensions)

    #Calculate local buckling at the TopHat punchout.
    top_hat_purlin_line.local_buckling_xx_net_pos = calculate_net_section_local_buckling_properties(top_hat_purlin_line)

    #Calculate distortional buckling including the influence of the TopHat punchout.
    top_hat_purlin_line.top_hat_purlin_distortional_net_cross_section_data = define_top_hat_purlin_distortional_net_section(top_hat_purlin_line)

    top_hat_purlin_line.distortional_buckling_xx_net_pos = calculate_net_section_distortional_buckling_properties(top_hat_purlin_line)

    #Calculate the first yield flexural strengths for each purlin line segment.  
    top_hat_purlin_line.yielding_flexural_strength_xx, top_hat_purlin_line.yielding_flexural_strength_xx_net, top_hat_purlin_line.yielding_flexural_strength_yy, top_hat_purlin_line.yielding_flexural_strength_free_flange_yy = calculate_yielding_flexural_strength(top_hat_purlin_line)

    #Calculate the local-global flexural strengths for each purlin line segment.   
    top_hat_purlin_line.local_global_flexural_strength_xx_no_hole, top_hat_purlin_line.local_global_flexural_strength_xx_hole, top_hat_purlin_line.local_global_flexural_strength_xx, top_hat_purlin_line.local_global_flexural_strength_yy, top_hat_purlin_line.local_global_flexural_strength_free_flange_yy = calculate_local_global_flexural_strength(top_hat_purlin_line)

    #Calculate distortional buckling strengths for each purlin line segment.
    top_hat_purlin_line.distortional_flexural_strength_xx = calculate_distortional_flexural_strength(top_hat_purlin_line)

    #Calculate torsion strength for each purlin line segment.
    top_hat_purlin_line.torsion_strength = calculate_torsion_strength(top_hat_purlin_line)

    #Calculate shear strength for each purlin line segment.
    top_hat_purlin_line.shear_strength_purlin, top_hat_purlin_line.shear_strength_top_hat, top_hat_purlin_line.shear_strength = calculate_shear_strength(top_hat_purlin_line)

    #Calculate web crippling strength at each support.
    #Assume purlin is limiting member here.  Don't consider TopHat.
    top_hat_purlin_line.web_crippling = PurlinLine.calculate_web_crippling_strength(purlin_line)


    return top_hat_purlin_line

end


function thin_walled_beam_interface(top_hat_purlin_line)

    #Discretize purlin line.
    member_definitions, dz, z, m = PurlinLine.discretize_purlin_line(top_hat_purlin_line)


    #Define ThinWalledBeam section property inputs.
    #Ix Iy Ixy J Cw

    num_purlin_sections = size(top_hat_purlin_line.inputs.purlin_cross_section_dimensions)[1]
    section_properties = Vector{Tuple{Float64, Float64, Float64, Float64, Float64}}(undef, num_purlin_sections)

    for i = 1:num_purlin_sections

        #Note -Ixy here since +y in ThinWalledBeam formulation is pointing down.
        section_properties[i] = (top_hat_purlin_line.top_hat_purlin_cross_section_data[i].section_properties.Ixx, top_hat_purlin_line.top_hat_purlin_cross_section_data[i].section_properties.Iyy, -top_hat_purlin_line.top_hat_purlin_cross_section_data[i].section_properties.Ixy, top_hat_purlin_line.top_hat_purlin_cross_section_data[i].section_properties.J, top_hat_purlin_line.top_hat_purlin_cross_section_data[i].section_properties.Cw)

    end


    #Define ThinWalledBeam material property inputs.
    num_purlin_materials = size(top_hat_purlin_line.inputs.purlin_material_properties)[1]

    material_properties = Vector{Tuple{Float64, Float64}}(undef, num_purlin_materials)

    for i = 1:num_purlin_materials

        material_properties[i] = (top_hat_purlin_line.inputs.purlin_material_properties[i][1], top_hat_purlin_line.inputs.purlin_material_properties[i][2])

    end

    #Define the lateral and rotational stiffness magnitudes for ThinWalledBeam.


    #There will be two sets of springs, one for the existing deck and one for the new deck.

    kx = Array{Array{Float64}}(undef, 2)
    kϕ = Array{Array{Float64}}(undef, 2)

    #First, the existing deck.

    num_purlin_segments = size(top_hat_purlin_line.bracing_data)[1]

    kx_segments = Vector{Float64}(undef, num_purlin_segments)
    kϕ_segments = Vector{Float64}(undef, num_purlin_segments)

    for i=1:num_purlin_segments

        kx_segments[i] = top_hat_purlin_line.bracing_data[i].kx 
        kϕ_segments[i]  = top_hat_purlin_line.bracing_data[i].kϕ
 
    end

    num_nodes = length(z)
    kx[1] = zeros(Float64, num_nodes)
    kϕ[1] = zeros(Float64, num_nodes)
    kx[1] .= kx_segments[m]
    kϕ[1] .= kϕ_segments[m]

    #Define the lateral spring location for ThinWalledBeam.    

    #Calculate the y-distance from the TopHat + purlin shear center to the existing deck lateral translational spring.
    spring_location_segment = Vector{Float64}(undef, num_purlin_sections)

    for i = 1:num_purlin_sections

        ys = top_hat_purlin_line.top_hat_purlin_cross_section_data[i].section_properties.ys  #distance from bottom fiber of purlin to shear center
        h = top_hat_purlin_line.inputs.purlin_cross_section_dimensions[i][5] 

        spring_location_segment[i] = h - ys  

    end


    ay_kx = Array{Array{Float64}}(undef, 2)

    #Define location of translational spring at each node.
    ay_kx[1] = Mesh.create_line_element_property_array(member_definitions, m, dz, spring_location_segment, 3, 1)


    #Now work on the new deck springs and spring locations.


    kx_segments = Vector{Float64}(undef, num_purlin_segments)
    kϕ_segments = Vector{Float64}(undef, num_purlin_segments)

    for i=1:num_purlin_segments

        kx_segments[i] = top_hat_purlin_line.new_deck_bracing_data[i].kx 
        kϕ_segments[i]  = top_hat_purlin_line.new_deck_bracing_data[i].kϕ
 
    end

    num_nodes = length(z)
    kx[2] = zeros(Float64, num_nodes)
    kϕ[2] = zeros(Float64, num_nodes)
    kx[2] .= kx_segments[m]
    kϕ[2] .= kϕ_segments[m]

    #Define the lateral spring location for ThinWalledBeam.    

    #Calculate the y-distance from the TopHat + purlin shear center to the new deck lateral translational spring.
    spring_location_segment = Vector{Float64}(undef, num_purlin_sections)

    for i = 1:num_purlin_sections

        ys = top_hat_purlin_line.top_hat_purlin_cross_section_data[i].section_properties.ys  #distance from bottom fiber of purlin to shear center
        top_hat_purlin_depth = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[i][4] + top_hat_purlin_line.inputs.purlin_cross_section_dimensions[i][5]
        #out to out TopHat + purlin height

        spring_location_segment[i] = top_hat_purlin_depth - ys  

    end

    #Define location of translational spring at each node.
    ay_kx[2] = Mesh.create_line_element_property_array(member_definitions, m, dz, spring_location_segment, 3, 1)


    #Define purlin line support locations for ThinWalledBeam.
    #location where u=v=ϕ=0
    supports = top_hat_purlin_line.inputs.support_locations

    #Define purlin line end boundary conditions for ThinWalledBeam.

    end_boundary_conditions = Array{Int64}(undef, 2)

    purlin_line_length = sum([top_hat_purlin_line.inputs.segments[i][1] for i=1:size(top_hat_purlin_line.inputs.segments)[1]])

    #type=1 u''=v''=ϕ''=0 (simply supported), type=2 u'=v'=ϕ'=0  (fixed), type=3 u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free end, e.g., a cantilever)

    #z=0 (left) end
    if supports[1] == 0.0
        end_boundary_conditions[1] = 1 #pin
    else
        end_boundary_conditions[1] = 3  #cantilever
    end

    #z=purlin_line_length (right) end
    if supports[end] == purlin_line_length
        end_boundary_conditions[2] = 1
    else
        end_boundary_conditions[2] = 3  #cantilever
    end

    #Calculate load magnitudes from user-defined pressure for ThinWalledBeam.
    q = top_hat_purlin_line.applied_pressure * top_hat_purlin_line.inputs.spacing #go from pressure to line load

    num_nodes = length(z)

    if q<0   #uplift wind pressure
        qx = zeros(num_nodes)
        qy = q .* ones(Float64, num_nodes)
    elseif q>= 0 #gravity pressure
        qx = -q .* sin(deg2rad(top_hat_purlin_line.inputs.roof_slope)) .* ones(Float64, num_nodes)
        qy = q .* cos(deg2rad(top_hat_purlin_line.inputs.roof_slope)) .* ones(Float64, num_nodes)
    end

    #Calculate the load locations for ThinWalledBeam, from the TopHat + purlin shear center.  
    ax_purlin_section = Vector{Float64}(undef, num_purlin_sections)
    ay_purlin_section = Vector{Float64}(undef, num_purlin_sections)

    for i = 1:num_purlin_sections

        center_top_flange_node_index = sum(top_hat_purlin_line.purlin_cross_section_data[i].n[1:3]) + sum(top_hat_purlin_line.purlin_cross_section_data[i].n_radius[1:3]) + floor(Int,top_hat_purlin_line.purlin_cross_section_data[i].n[4]/2) + 1

        ax_purlin_section[i] = top_hat_purlin_line.purlin_cross_section_data[i].node_geometry[center_top_flange_node_index, 1] - top_hat_purlin_line.top_hat_purlin_cross_section_data[i].section_properties.xs

        top_hat_purlin_depth = top_hat_purlin_line.inputs.top_hat_cross_section_dimensions[i][4] + top_hat_purlin_line.inputs.purlin_cross_section_dimensions[i][5]

        #This is different since the load is now applied at the top of the TopHat.
        ay_purlin_section[i] = top_hat_purlin_depth - top_hat_purlin_line.top_hat_purlin_cross_section_data[i].section_properties.ys
        
    end

    #Define the load location at each node.
    ax = Mesh.create_line_element_property_array(member_definitions, m, dz, ax_purlin_section, 3, 1)
    ay = Mesh.create_line_element_property_array(member_definitions, m, dz, ay_purlin_section, 3, 1)

    return z, m, member_definitions, section_properties, material_properties, kx, kϕ, ay_kx, qx, qy, ax, ay, end_boundary_conditions, supports

end



function beam_column_interface(top_hat_purlin_line)

    #Discretize purlin line.
    member_definitions, dz, z, m = PurlinLine.discretize_purlin_line(top_hat_purlin_line)

    #Define the number of purlin cross-sections.
    num_purlin_sections = size(top_hat_purlin_line.inputs.purlin_cross_section_dimensions)[1]

    #Initialize an array of tuples to hold the free flange section properties.
    section_properties = Vector{Tuple{Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64,}}(undef, num_purlin_sections)

    for i = 1:num_purlin_sections

        Af = top_hat_purlin_line.free_flange_cross_section_data[i].section_properties.A
        Ixf = top_hat_purlin_line.free_flange_cross_section_data[i].section_properties.Ixx
        Iyf = top_hat_purlin_line.free_flange_cross_section_data[i].section_properties.Iyy
        Jf = top_hat_purlin_line.free_flange_cross_section_data[i].section_properties.J
        Cwf = top_hat_purlin_line.free_flange_cross_section_data[i].section_properties.Cw
        xcf = top_hat_purlin_line.free_flange_cross_section_data[i].section_properties.xc
        ycf = top_hat_purlin_line.free_flange_cross_section_data[i].section_properties.yc
        xsf = top_hat_purlin_line.free_flange_cross_section_data[i].section_properties.xs
        ysf = top_hat_purlin_line.free_flange_cross_section_data[i].section_properties.ys

        section_properties[i] = (Af, Ixf, Iyf, Jf, Cwf, xcf, ycf, xsf, ysf)

    end


    #Define BeamColumn material property inputs.
    num_purlin_materials = size(top_hat_purlin_line.inputs.purlin_material_properties)[1]

    material_properties = Vector{Tuple{Float64, Float64}}(undef, num_purlin_materials)

    for i = 1:num_purlin_materials

        material_properties[i] = (top_hat_purlin_line.inputs.purlin_material_properties[i][1], top_hat_purlin_line.inputs.purlin_material_properties[i][2])

    end

   
    num_purlin_segments = size(top_hat_purlin_line.bracing_data)[1]

    #Define kxf along the purlin line.
    kxf_segments = [top_hat_purlin_line.free_flange_data[i].kxf for i=1:num_purlin_segments]

    num_nodes = length(z)
    kxf = zeros(Float64, num_nodes)
    kxf .= kxf_segments[m]

    #There is no kyf assumed.
    kyf = zeros(Float64, num_nodes)

    #Define kϕf along the purlin line.
    kϕf_segments = [top_hat_purlin_line.free_flange_data[i].kϕf for i=1:num_purlin_segments]
    kϕf = zeros(Float64, num_nodes)
    kϕf .= kϕf_segments[m]    

    #Assume the lateral spring acts at the free flange centroid.  This means hx =hy = 0.
    hxf = zeros(Float64, num_nodes)
    hyf = zeros(Float64, num_nodes)

    #Define shear flow force in free flange.

    #Define the purlin segment properties.
    kH_segments = [top_hat_purlin_line.free_flange_data[i].kH for i=1:num_purlin_segments]
    kH = zeros(Float64, num_nodes)
    kH .= kH_segments[m]

    #The shear flow is applied at the free flange centerline.  The distance ay in StructuresKit.BeamColumn is the distance from the shear center to the load along the centroidal y-axis.   Since the shear center for just the free flange is close to the free flange centerline, assume ay= 0.  

    ayf = zeros(Float64, num_nodes)

    #There is no qyf so this can be set to zero.
    axf = zeros(Float64, num_nodes)

    #Define supports.   Combine frame supports and intermediate bridging here.
    supports = sort(unique([top_hat_purlin_line.inputs.support_locations; top_hat_purlin_line.inputs.bridging_locations]))

    #Define purlin line end boundary conditions for BeamColumn.

    end_boundary_conditions = Array{Int64}(undef, 2)

    purlin_line_length = sum([top_hat_purlin_line.inputs.segments[i][1] for i=1:size(top_hat_purlin_line.inputs.segments)[1]])

    #type=1 u''=v''=ϕ''=0 (simply supported), type=2 u'=v'=ϕ'=0  (fixed), type=3 u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free end, e.g., a cantilever)

    #z=0 (left) end
    if supports[1] == 0.0
        end_boundary_conditions[1] = 1 #pin
    else
        end_boundary_conditions[1] = 3  #cantilever
    end

    #z=purlin_line_length (right) end
    if supports[end] == purlin_line_length
        end_boundary_conditions[2] = 1
    else
        end_boundary_conditions[2] = 3  #cantilever
    end


    return z, m, member_definitions, section_properties, material_properties, kxf, kyf, kϕf, kH, hxf, hyf, axf, ayf, end_boundary_conditions, supports

end


function calculate_free_flange_axial_force(Mxx, member_definitions, top_hat_purlin_line)
    #this may not need to be updated, could use PurlinLine version
    num_purlin_sections = size(top_hat_purlin_line.inputs.purlin_cross_section_dimensions)[1]

    P_unit = zeros(Float64, num_purlin_sections)

    #Loop over the purlin cross-sections in the line.
    for i = 1:num_purlin_sections

        #Find web node at H/5.
        web_index = top_hat_purlin_line.purlin_cross_section_data[i].n[1] + top_hat_purlin_line.purlin_cross_section_data[i].n_radius[1] + top_hat_purlin_line.purlin_cross_section_data[i].n[2] + top_hat_purlin_line.purlin_cross_section_data[i].n_radius[2] + 1 + 1

        #Use the local_buckling_xx_pos node geometry and reference stress from CUFSM (Mxx = 1).
        dx = diff(top_hat_purlin_line.local_buckling_xx_pos[i].CUFSM_data.node[1:web_index,2])
        dy = diff(top_hat_purlin_line.local_buckling_xx_pos[i].CUFSM_data.node[1:web_index,3])
        ds = sqrt.(dx.^2 .+ dy.^2)
        s = [0; cumsum(ds)]   #line coordinates around free flange

        #Integrate the reference stress (Mxx = 1.0) in the free flange to find the reference axial force.   
        stress = top_hat_purlin_line.local_buckling_xx_pos[i].CUFSM_data.node[1:web_index,8] 
        t = top_hat_purlin_line.inputs.purlin_cross_section_dimensions[i][2]
        P_unit[i] = integrate(s, stress) * t

    end

    #Scale the reference axial force along the purlin line to define the axial force in the free flange.
    #The sign convention for P is + (compression), - (tension) to match StructuresKit.BeamColumn.
    dz = diff(top_hat_purlin_line.model.z)
    P = Mesh.create_line_element_property_array(member_definitions, top_hat_purlin_line.model.m, dz, P_unit, 3, 1) .* Mxx

    return P

end


function analysis(top_hat_purlin_line)

    z, m, member_definitions, section_properties, material_properties, kx, kϕ, ay_kx, qx, qy, ax, ay, end_boundary_conditions, supports = TopHatDesigner.thin_walled_beam_interface(top_hat_purlin_line)

    #Set up ThinWalledBeam model.
    model = ThinWalledBeam.define(z, m, member_definitions, section_properties, material_properties, kx, kϕ, ay_kx, qx, qy, ax, ay, end_boundary_conditions, supports)

    #Solve ThinWalledBeam model.
    top_hat_purlin_line.model = ThinWalledBeam.solve(model)

    #Calculate purlin line internal forces and moments from deformations, add them to data structure.
    Mxx, Myy, Vxx, Vyy, T, B = PurlinLine.calculate_internal_forces(top_hat_purlin_line.model)

    num_nodes = length(top_hat_purlin_line.model.z)
    P = zeros(Float64, num_nodes)  #No axial force in purlin for now.  Could be added later.

    #Add internal forces to data structure.
    top_hat_purlin_line.internal_forces = PurlinLine.InternalForceData(P, Mxx, Myy, Vxx, Vyy, T, B)

    #Translate purlin_line design variables to BeamColumn design variables.
    z, m, member_definitions, section_properties, material_properties, kxf, kyf, kϕf, kH, hxf, hyf, axf, ayf, end_boundary_conditions, supports = beam_column_interface(top_hat_purlin_line)

    #Calculate axial force in free flange.
    Pf = calculate_free_flange_axial_force(Mxx, member_definitions, top_hat_purlin_line)

    #Apply the shear flow based on the y-direction load along the purlin line free flange model.
    qxf = qy .* kH

    #The y-direction load is assumed to be zero in the free flange model.
    num_nodes = length(z)
    qyf = zeros(Float64, num_nodes)

    #Set up the free flange model.
    top_hat_purlin_line.free_flange_model = BeamColumn.define(z, m, member_definitions, section_properties, material_properties, kxf, kyf, kϕf, hxf, hyf, qxf, qyf, Pf, axf, ayf, end_boundary_conditions, supports)

    #Run the free flange model.
    top_hat_purlin_line.free_flange_model = BeamColumn.solve(top_hat_purlin_line.free_flange_model)

    #Calculate internal forces in the free flange.
    Mxx, Myy, Vxx, Vyy, T, B = PurlinLine.calculate_internal_forces(top_hat_purlin_line.free_flange_model)

    #Add free flange internal forces to data structure.
    top_hat_purlin_line.free_flange_internal_forces = PurlinLine.InternalForceData(Pf, Mxx, Myy, Vxx, Vyy, T, B)



    #Calculate demand-to-capacity ratios for each of the purlin line limit states.
    top_hat_purlin_line.flexure_torsion_demand_to_capacity, eMnℓ_xx_all, eMnℓ_yy_all, eBn_all, eMnℓ_yy_free_flange_all = PurlinLine.calculate_flexure_torsion_demand_to_capacity(top_hat_purlin_line)
    top_hat_purlin_line.distortional_demand_to_capacity, eMnd_xx_all = PurlinLine.calculate_distortional_buckling_demand_to_capacity(top_hat_purlin_line)
    top_hat_purlin_line.flexure_shear_demand_to_capacity, eMnℓ_xx_all, eVn_all = PurlinLine.calculate_flexure_shear_demand_to_capacity(top_hat_purlin_line)        
    top_hat_purlin_line.biaxial_bending_demand_to_capacity, eMnℓ_xx_all, eMnℓ_yy_all = PurlinLine.calculate_biaxial_bending_demand_to_capacity(top_hat_purlin_line)
    top_hat_purlin_line.web_crippling_demand_to_capacity = PurlinLine.calculate_web_crippling_demand_to_capacity(top_hat_purlin_line)

    #Add expected strengths along purlin line to data structure.
    top_hat_purlin_line.expected_strengths = PurlinLine.ExpectedStrengths(eMnℓ_xx_all, eMnℓ_yy_all, eMnℓ_yy_free_flange_all, eMnd_xx_all, eVn_all, eBn_all)

    return top_hat_purlin_line

end

function capacity(top_hat_purlin_line)

    DC_tolerance = 0.01  
    
    if top_hat_purlin_line.loading_direction == "gravity"

        load_sign = 1.0
    
    elseif top_hat_purlin_line.loading_direction == "uplift"
    
        load_sign = -1.0
    
    end

    #Run a very small pressure to get the test going.
    top_hat_purlin_line.applied_pressure = load_sign * 10^-6
    top_hat_purlin_line = TopHatDesigner.analysis(top_hat_purlin_line)
    max_DC = PurlinLine.find_max_demand_to_capacity(top_hat_purlin_line)

    #Define initial residual.
    residual = 1.0 - abs(max_DC)

    while residual > DC_tolerance

        new_pressure = top_hat_purlin_line.applied_pressure / max_DC
        top_hat_purlin_line.applied_pressure = top_hat_purlin_line.applied_pressure + (new_pressure - top_hat_purlin_line.applied_pressure) / 2

        top_hat_purlin_line = TopHatDesigner.analysis(top_hat_purlin_line)
        max_DC = PurlinLine.find_max_demand_to_capacity(top_hat_purlin_line)

        residual = 1.0 - abs(max_DC)

    end

    top_hat_purlin_line.failure_limit_state, top_hat_purlin_line.failure_location = PurlinLine.identify_failure_limit_state(top_hat_purlin_line)


    return top_hat_purlin_line

end

end # module
