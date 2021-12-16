# ENGRI-1120-Project
### A Pluto.jl notebook ###
# v0.17.2
 
using Markdown
using InteractiveUtils
 
# ╔═╡ 5458cafc-430a-4e2e-a3f9-d23023e6053b
begin
 
        # load some external packages 
        using PlutoUI
        using DataFrames
        using BSON
        using GLPK
        using PrettyTables
        using Plots
        
        # setup my paths (where are my files?)
        _PATH_TO_ROOT = pwd() 
        _PATH_TO_SRC = joinpath(_PATH_TO_ROOT,"src")
        _PATH_TO_MODEL = joinpath(_PATH_TO_ROOT,"model")
        _PATH_TO_FIGS = joinpath(_PATH_TO_ROOT,"figs")
        
        # load the ENGRI 1120 project code library -
        include(joinpath(_PATH_TO_SRC,"Include.jl"))
 
        # load the model -
        MODEL = BSON.load(joinpath(_PATH_TO_MODEL,"model_v2.bson"), @__MODULE__)
 
        # show -
        nothing
end
 
# ╔═╡ 6f9851d0-9573-4965-8456-92d2c57afa91
md"""
## ENGRI 1120: Design and Analysis of a Sustainable Cell-Free Production Process for Industrially Important Small Molecules
"""
 
# ╔═╡ 833d7250-f4ab-49a2-a43e-1b21def59ad4
html"""
<p style="font-size:20px;">Creative Name: Sebastian Soto, Joan Varous, Ashlyn Dumaw </br>
Smith School of Chemical and Biomolecular Engineering, Cornell University, Ithaca NY 14850</p>
"""
 
# ╔═╡ 251363ad-1927-4b05-99f5-8c3f2508c0cb
md"""
### Introduction
"""
 
# ╔═╡ 8d20604d-360a-4524-bf38-1ac30e5d9be4
md"""
For our project we had to maximize the output of PDGN using the sugar Xylose. We want to see if this process would be an effective and profitable process.
"""
 
# ╔═╡ 884a0a7d-e5d8-4417-b109-d00c37a01766
md"""
### Materials and Methods
"""
 
# ╔═╡ e9ec4f6d-aea4-4ff4-b25e-49f422a4051f
md"""
For our materials, we used D-Xylose, Potassium Nitrate, NADH, and Beta-D-Glucose 6-Phosphate. We also ended up needing 11 different chips and 110 separators. We also needed 11 pumps to pump in the material since we chose to run our reactors in parallel.
 
We maximized the output of PDGN using a parallel system. We also decided to separate out other profitable products by separating the outflow stream. When 100% of the stream was going into the separator, we still needed 11 reactors running simultaneously in order to produce a little over 1 g/hr of PDGN. Also, to get at least 95% purity, we needed at least 6 levels of separation. So, using this information, we found the minimum fraction of flow that needs to go into the first separator in order to still get at least 1 g/hr of PDGN after 6 levels of separation. That number turned out to be 92% of the outflow stream. Then after calculating that, we wrote another code to calculate the output of a different material from the 8% of the outflow stream remaining. The material that was produced that would’ve yielded us the most profit was NAD and it only required 4 levels of separation. On top of that, we also we calculated it so that we would use the minimum materials required to get the maximum output of PDGN possible. So all of our inputs get completely used up so they aren't present in the outflow stream that gets separated.
 
The reaction process we calculated should look like this:
 
1)D-Xylose + NADPH + H+ ⇔ Xylitol + NADP+ 
 
2)beta-D-Glucose 6-phosphate + NADP+ ⇔ D-Glucono-15-lactone 6-phosphate + NADPH + H+ 
 
3)Xylitol + NAD+ ⇔ D-Xylulose + NADH + H+ 
 
4)ATP + D-Xylulose ⇔ ADP + D-Xylulose 5-phosphate 
 
5)D-Xylulose 5-phosphate ⇔ D-Ribulose 5-phosphate 
 
6)D-Ribulose 5-phosphate ⇔ D-Ribose 5-phosphate 
 
7)D-Ribose 5-phosphate + D-Xylulose 5-phosphate ⇔ Sedoheptulose 7-phosphate + D-Glyceraldehyde 3-phosphate 
 
8)D-Glyceraldehyde 3-phosphate ⇔ Glycerone phosphate 
 
9)Glycerone phosphate ⇔ Methylglyoxal + Orthophosphate 
 
10)Methylglyoxal + NADPH + H+ ⇔ (S)-Lactaldehyde + NADP+ 
 
11)Methylglyoxal + NADH + H+ ⇔ (S)-Lactaldehyde + NAD 
 
12)(S)-Lactaldehyde + NADH + H+ ⇔ Propane-1,2-diol + NAD+ 
 
13)Propane-1,2-diol + 2KNO3 ⇔ PGDN + 2KOH 
 
"""
 
# ╔═╡ ad5d595e-4dba-49cd-a446-e1df737fd75d
md"""
##### Step 1: Configure the Flux Balance Analysis (FBA) calculation for a _single_ chip
"""
 
# ╔═╡ 5bfdf6f9-2927-4a9a-a386-8840c676329b
begin
 
        # setup the FBA calculation for the project -
 
        # === SELECT YOUR PRODUCT HERE ==================================================== #
        # What rate are trying to maximize? (PDGN)
        # rn:R08199 = isoprene
        # rn:28235c0c-ec00-4a11-8acb-510b0f2e2687 = PGDN
        # rn:rn:R09799 = Hydrazine
        # rn:R03119 = 3G
        idx_target_rate = find_reaction_index(MODEL,:reaction_number=>"rn:28235c0c-ec00-4a11-8acb-510b0f2e2687")
        # ================================================================================= #
 
        # First, let's build the stoichiometric matrix from the model object -
        (cia,ria,S) = build_stoichiometric_matrix(MODEL);
 
        # Next, what is the size of the system? (ℳ = number of metabolites, ℛ = number of reactions)
        (ℳ,ℛ) = size(S)
 
        # Next, setup a default bounds array => update specific elements
        # We'll correct the directionality below -
        Vₘ = (13.7)*(3600)*(50e-9)*(1000) # units: mmol/hr
        flux_bounds = [-Vₘ*ones(ℛ,1) Vₘ*ones(ℛ,1)]
 
        # update the flux bounds -> which fluxes can can backwards? 
        # do determine this: sgn(v) = -1*sgn(ΔG)
        updated_flux_bounds = update_flux_bounds_directionality(MODEL,flux_bounds)
 
        # hard code some bounds that we know -
        updated_flux_bounds[44,1] = 0.0  # ATP synthesis can't run backwards 
 
        # What is the default mol flow input array => update specific elements
        # strategy: start with nothing in both streams, add material(s) back
        n_dot_input_stream_1 = zeros(ℳ,1)        # stream 1
        n_dot_input_stream_2 = zeros(ℳ,1)        # stream 2
 
        # === YOU NEED TO CHANGE BELOW HERE ====================================================== #
        # Let's lookup stuff that we want/need to supply to the chip to get the reactiont to go -
        # what you feed *depends upon your product*
        compounds_that_we_need_to_supply_feed_1 = [
               "d-xylose", 
        ]
 
        # what are the amounts that we need to supply to chip in feed stream 1 (units: mmol/hr)?
        mol_flow_values_feed_1 = [
               0.634         ; # d-xylose mmol/hr
                
        ]
 
        # what is coming into feed stream 2?
        compounds_that_we_need_to_supply_feed_2 = [
                "nadh", "potassium nitrate", "beta-d-glucose 6-phosphate"
        ]
 
        # let's always add Vₘ into feed stream 2
        mol_flow_values_feed_2 = [
                2.11408                ; # nadh mmol/hr
               4.932 ; # potassium nitrate mmol/hr
                1.05686 ; # beta-d-glucose 6-phosphate mmol/hr
        ]
        
        
        # === YOU NEED TO CHANGE ABOVE HERE ====================================================== #
 
        # stream 1:
        idx_supply_stream_1 = Array{Int64,1}()
        for compound in compounds_that_we_need_to_supply_feed_1
               idx = find_compound_index(MODEL,:compound_name=>compound)
                push!(idx_supply_stream_1,idx)
        end
 
        # stream 2:
        idx_supply_stream_2 = Array{Int64,1}()
        for compound in compounds_that_we_need_to_supply_feed_2
               idx = find_compound_index(MODEL,:compound_name=>compound)
                push!(idx_supply_stream_2,idx)
        end
        
        # supply for stream 1 and stream 2
        n_dot_input_stream_1[idx_supply_stream_1] .= mol_flow_values_feed_1
        n_dot_input_stream_2[idx_supply_stream_2] .= mol_flow_values_feed_2
        
        # setup the species bounds array -
        species_bounds = [-1.0*(n_dot_input_stream_1.+n_dot_input_stream_2) 1000.0*ones(ℳ,1)]
 
        # Lastly, let's setup the objective function -
        c = zeros(ℛ)
        c[idx_target_rate] = -1.0
 
        # show -
        nothing
end
 
# ╔═╡ e8a4faf8-2285-4544-830c-f39d3847e8cc
md"""
##### Step 2: Method to compute the composition that is going into the downstream separation system 
 
In a parallel system of $N$ chips, each chip acts independently. Thus, to compute the output from the mixer operation we add up the components in each stream into the mixer (N inputs and a single output). Starting from the steady-state species mol balance: 
 
$$\sum_{s=1}^{N+1}v_{s}\dot{n}_{i,s} = 0\qquad{i=1,2,\dots,\mathcal{M}}$$
 
we can solve for the mixer output stream composition:
 
$$\dot{n}_{i,N+1} = -\sum_{s=1}^{N}v_{s}\dot{n}_{i,s}\qquad{i=1,2,\dots,\mathcal{M}}$$
 
However, since each chip is _identical_ we know that: $\dot{n}_{i,N+1} = N\times\dot{n}_{i,1}$. Alternatively, we could do the same thing with species mass balances (which is probably more useful in this case because our downstream separation units operate on a mass basis).
 
"""
 
# ╔═╡ 10424555-39cc-4ddf-8c22-db91cf102bfd
md"""
##### Step 3: Downstream separation using Magical Sepration Units (MSUs)
 
To separate the desired product from the unreacted starting materials and by-products, let's suppose the teaching team invented a magical separation unit or MSU. MSUs have one stream in, and two streams out (called the top, and bottom, respectively) and a fixed separation ratio for all products (that's what makes them magical), where the desired product is _always_ in the top stream at some ratio $\theta$. In particular, if we denote $i=\star$ as the index for the desired product (in this case 1,3 propanediol), then after one pass (stream 1 is the input, stream 2 is the top, and stream 3 is the bottom) we have:
 
$$\begin{eqnarray}
\dot{m}_{\star,2} &=& \theta_{\star}\dot{m}_{\star,1}\\
\dot{m}_{\star,3} &=& (1-\theta_{\star})\dot{m}_{\star,1}\\
\end{eqnarray}$$
 
for the product. In this case, we set $\theta_{\star}$ = 0.75. On the other hand, for _all_ other materials in the input, we have $\left(1-\theta_{\star}\right)$ in the top, and $\theta_{\star}$ in the bottom, i.e.,
 
$$\begin{eqnarray}
\dot{m}_{i,2} &=& (1-\theta_{\star})\dot{m}_{i,1}\qquad{\forall{i}\neq\star}\\
\dot{m}_{i,3} &=& \theta_{\star}\dot{m}_{i,1}\\
\end{eqnarray}$$
 
If we chain these units together we can achieve a desired degree of separation.
"""
 
# ╔═╡ a3632011-833e-431b-b08f-f2896ad0a82a
md"""
### Results and Discussion
"""
 
# ╔═╡ 39f3633e-b9df-4d01-946f-ba2d6c8ba6b7
begin
 
        # compute the optimal flux -
        result = calculate_optimal_flux_distribution(S, updated_flux_bounds, species_bounds, c);
 
        # get the open extent vector -
        ϵ_dot = result.calculated_flux_array
 
        # what is the composition coming out of the first chip?
        n_dot_out_chip_1 = (n_dot_input_stream_1 + n_dot_input_stream_2 + S*ϵ_dot);
 
        # did this converge?
        with_terminal() do
 
               # get exit/status information from the solver -
                exit_flag = result.exit_flag
                status_flag = result.status_flag
 
               # display -
                println("Computed optimal flux distribution w/exit_flag = 0: $(exit_flag==0) and status_flag = 5: $(status_flag == 5)")
        end
end
 
# ╔═╡ 4b3ef98c-d304-4ef4-95ef-f1d1ce562e36
md"""
###### Table 1: State table from a single chip (species mol flow rate mmol/hr at exit)
"""
 
# ╔═╡ 7166a917-b676-465c-a441-4ff0530faf92
begin
 
        # compute the mol flow rate out of the device -
        n_dot_output = (n_dot_input_stream_1 + n_dot_input_stream_2 + S*ϵ_dot);
 
        # get the array of MW -
        MW_array = MODEL[:compounds][!,:compound_mw]
 
        # convert the output mol stream to a mass stream -
        mass_dot_output = (n_dot_output.*MW_array)*(1/1000)
 
        # what is the total coming out?
        total_mass_out = sum(mass_dot_output)
        
        # display -
        with_terminal() do
 
               # what are the compound names and code strings? -> we can get these from the MODEL object 
                compound_name_strings = MODEL[:compounds][!,:compound_name]
                compound_id_strings = MODEL[:compounds][!,:compound_id]
               
                # how many molecules are in the state array?
                ℳ_local = length(compound_id_strings)
        
               # initialize some storage -
                state_table = Array{Any,2}(undef,ℳ_local,9)
 
               # get the uptake array from the result -
                uptake_array = result.uptake_array
 
               # populate the state table -
               for compound_index = 1:ℳ_local
                       state_table[compound_index,1] = compound_index
                       state_table[compound_index,2] = compound_name_strings[compound_index]
                       state_table[compound_index,3] = compound_id_strings[compound_index]
                       state_table[compound_index,4] = n_dot_input_stream_1[compound_index]
                       state_table[compound_index,5] = n_dot_input_stream_2[compound_index]
                       
 
                       # for display -
                       tmp_value = abs(n_dot_output[compound_index])
                        state_table[compound_index,6] = (tmp_value) <= 1e-6 ? 0.0 : n_dot_output[compound_index]
 
                       # show the Δ -
                       tmp_value = abs(uptake_array[compound_index])
                       state_table[compound_index,7] = (tmp_value) <= 1e-6 ? 0.0 : uptake_array[compound_index]
 
                       # show the mass -
                       tmp_value = abs(mass_dot_output[compound_index])
                       state_table[compound_index,8] = (tmp_value) <= 1e-6 ? 0.0 : mass_dot_output[compound_index]
 
                       # show the mass fraction -
                       # show the mass -
                       tmp_value = abs(mass_dot_output[compound_index])
                       state_table[compound_index,9] = (tmp_value) <= 1e-6 ? 0.0 : (1/total_mass_out)*mass_dot_output[compound_index]
               end
 
               # header row -
                state_table_header_row = (["i","name","id","n₁_dot", "n₂_dot", "n₃_dot","Δn_dot", "m₃_dot", "ωᵢ_output"],
                       ["","","","mmol/hr", "mmol/hr", "mmol/hr", "mmol/hr", "g/hr", ""]);
               
               # write the table -
                pretty_table(state_table; header=state_table_header_row)
        end
end
 
# ╔═╡ 80205dc2-0cd9-4543-be6c-2b3a7a5010d5
# How many chips do we want to operate in parallel?
number_of_chips = 11;
 
# ╔═╡ eb091c37-29f6-45e8-8716-126c2df7f125
# how many levels are we going to have in the separation tree?
number_of_levels = 6;
 
# ╔═╡ 65c26314-f7de-42c7-978c-5fe18ef45850
# what compound are we trying to separate?
idx_target_compound = find_compound_index(MODEL,:compound_name=>"PGDN");
 
# ╔═╡ fe1a84e2-0a44-4341-9add-35f8bb296454
begin
 
        # alias the mass flow into the sep-units
        # mass flow coming out of the mixer -
        mass_flow_into_seps = (number_of_chips)*(.92)mass_dot_output
        
        # define the split -
        θ = 0.75
 
        # most of the "stuff" has a 1 - θ in the up, and a θ in the down
        u = (1-θ)*ones(ℳ,1)
        d = θ*ones(ℳ,1)
 
        # correct defaults -
        u[idx_target_compound] = θ
        d[idx_target_compound] = 1 - θ
 
        # let's compute the composition of the *always up* stream -
        
        # initialize some storage -
        species_mass_flow_array_top = zeros(ℳ,number_of_levels)
        species_mass_flow_array_bottom = zeros(ℳ,number_of_levels)
 
        for species_index = 1:ℳ
               value = mass_flow_into_seps[species_index]
                species_mass_flow_array_top[species_index,1] = value
                species_mass_flow_array_bottom[species_index,1] = value
        end
        
        for level = 2:number_of_levels
 
               # compute the mass flows coming out of the top -
                m_dot_top = mass_flow_into_seps.*(u.^(level-1))
                m_dot_bottom = mass_flow_into_seps.*(d.^(level-1))
 
               # update my storage array -
               for species_index = 1:ℳ
                       species_mass_flow_array_top[species_index,level] = m_dot_top[species_index]
                       species_mass_flow_array_bottom[species_index,level] = m_dot_bottom[species_index]
               end
        end
        
        # what is the mass fraction in the top stream -
        species_mass_fraction_array_top = zeros(ℳ,number_of_levels)
        species_mass_fraction_array_bottom = zeros(ℳ,number_of_levels)
 
        # array to hold the *total* mass flow rate -
        total_mdot_top_array = zeros(number_of_levels)
        total_mdot_bottom_array = zeros(number_of_levels)
        
        # this is a dumb way to do this .. you're better than that JV come on ...
        T_top = sum(species_mass_flow_array_top,dims=1)
        T_bottom = sum(species_mass_flow_array_bottom,dims=1)
        for level = 1:number_of_levels
 
               # get the total for this level -
                T_level_top = T_top[level]
                T_level_bottom = T_bottom[level]
 
               # grab -
                total_mdot_top_array[level] = T_level_top
                total_mdot_bottom_array[level] = T_level_bottom
 
               for species_index = 1:ℳ
                       species_mass_fraction_array_top[species_index,level] = (1/T_level_top)*
                               (species_mass_flow_array_top[species_index,level])
                       species_mass_fraction_array_bottom[species_index,level] = (1/T_level_bottom)*
                               (species_mass_flow_array_bottom[species_index,level])
               end
        end
end
 
# ╔═╡ efe968b6-4914-4c4c-a2fb-50d7e71f582b
begin
 
        stages = (1:number_of_levels) |> collect
        plot(stages,species_mass_fraction_array_top[idx_target_compound,:], linetype=:steppre,lw=2,legend=:bottomright, 
                label="Mass fraction i = PDO Tops")
        xlabel!("Stage index l",fontsize=18)
        ylabel!("Tops mass fraction ωᵢ (dimensionless)",fontsize=18)
 
        # make a 0.95 line target line -
        target_line = 0.95*ones(number_of_levels)
        plot!(stages, target_line, color="red", lw=2,linestyle=:dash, label="Target 95% purity")
end
 
# ╔═╡ 4a308e7a-0149-4816-a24e-ecf23c0a759c
with_terminal() do
 
        # initialize some space -
        state_table = Array{Any,2}(undef, number_of_levels, 3)
        for level_index = 1:number_of_levels
                state_table[level_index,1] = level_index
                state_table[level_index,2] = species_mass_fraction_array_top[idx_target_compound, level_index]
                state_table[level_index,3] = total_mdot_top_array[level_index]
        end
        
        # header -
        state_table_header_row = (["stage","ωᵢ i=⋆ top","mdot"],
                       ["","","g/hr"]);
 
        # write the table -
        pretty_table(state_table; header=state_table_header_row)
end
 
# ╔═╡ e4ae84ea-fa9e-4588-b986-fb779467ec88
md"""
### Seperating the output stream and using the MSU for a different species 
"""
 
# ╔═╡ f10812e0-1e33-44a6-88d9-bb252563e3bd
# what compound are we trying to separate?
idx_target_compound2 = find_compound_index(MODEL,:compound_name=>"nad");
 
# ╔═╡ 872c5395-047c-4c14-b37a-a68933879656
# how many levels are we going to have in the separation tree?
number_of_levels2 = 6;
 
# ╔═╡ bfb3f60d-f2bc-4f17-8407-f86291fd55b2
begin
 
        # alias the mass flow into the sep-units
        # mass flow coming out of the mixer -
        mass_flow_into_seps2 = (number_of_chips)*(.08)mass_dot_output
        
        # define the split -
        θ2 = 0.75
 
        # most of the "stuff" has a 1 - θ in the up, and a θ in the down
        u2 = (1-θ)*ones(ℳ,1)
        d2 = θ*ones(ℳ,1)
 
        # correct defaults -
        u2[idx_target_compound2] = θ2
        d2[idx_target_compound2] = 1 - θ2
 
        # let's compute the composition of the *always up* stream -
        
        # initialize some storage -
        species_mass_flow_array_top2 = zeros(ℳ,number_of_levels2)
        species_mass_flow_array_bottom2 = zeros(ℳ,number_of_levels2)
 
        for species_index = 1:ℳ
               value = mass_flow_into_seps[species_index]
                species_mass_flow_array_top2[species_index,1] = value
                species_mass_flow_array_bottom2[species_index,1] = value
        end
        
        for level = 2:number_of_levels
 
               # compute the mass flows coming out of the top -
                m_dot_top2 = mass_flow_into_seps2.*(u2.^(level-1))
                m_dot_bottom2 = mass_flow_into_seps2.*(d2.^(level-1))
 
               # update my storage array -
               for species_index = 1:ℳ
                       species_mass_flow_array_top2[species_index,level] = m_dot_top2[species_index]
                       species_mass_flow_array_bottom2[species_index,level] = m_dot_bottom2[species_index]
               end
        end
        
        # what is the mass fraction in the top stream -
        species_mass_fraction_array_top2 = zeros(ℳ,number_of_levels2)
        species_mass_fraction_array_bottom2 = zeros(ℳ,number_of_levels2)
 
        # array to hold the *total* mass flow rate -
        total_mdot_top_array2 = zeros(number_of_levels2)
        total_mdot_bottom_array2 = zeros(number_of_levels2)
        
        # this is a dumb way to do this ... you're better than that JV come on ...
        T_top2 = sum(species_mass_flow_array_top2,dims=1)
        T_bottom2 = sum(species_mass_flow_array_bottom2,dims=1)
        for level = 1:number_of_levels2
 
               # get the total for this level -
                T_level_top2 = T_top2[level]
                T_level_bottom2 = T_bottom2[level]
 
               # grab -
                total_mdot_top_array2[level] = T_level_top2
                total_mdot_bottom_array2[level] = T_level_bottom2
 
               for species_index = 1:ℳ
                       species_mass_fraction_array_top2[species_index,level] = (1/T_level_top2)*
                               (species_mass_flow_array_top2[species_index,level])
                       species_mass_fraction_array_bottom2[species_index,level] = (1/T_level_bottom2)*
                               (species_mass_flow_array_bottom2[species_index,level])
               end
        end
end
 
# ╔═╡ cd239ea5-3947-48e7-95e4-3db2c90a612f
begin
 
        stages2 = (1:number_of_levels2) |> collect
        plot(stages2,species_mass_fraction_array_top2[idx_target_compound2,:], linetype=:steppre,lw=2,legend=:bottomright, 
                label="Mass fraction i = PDO Tops")
        xlabel!("Stage index l",fontsize=18)
        ylabel!("Tops mass fraction ωᵢ (dimensionless)",fontsize=18)
 
        # make a 0.95 line target line -
        target_line2 = 0.95*ones(number_of_levels2)
        plot!(stages2, target_line2, color="red", lw=2,linestyle=:dash, label="Target 95% purity")
end
 
# ╔═╡ a349cd6b-f6d9-4544-a735-d503db4e818b
with_terminal() do
 
        # initialize some space -
        state_table = Array{Any,2}(undef, number_of_levels2, 3)
        for level_index = 1:number_of_levels2
                state_table[level_index,1] = level_index
                state_table[level_index,2] = species_mass_fraction_array_top2[idx_target_compound2, level_index]
                state_table[level_index,3] = total_mdot_top_array2[level_index]
        end
        
        # header -
        state_table_header_row = (["stage","ωᵢ i=⋆ top","mdot"],
                       ["","","g/hr"]);
 
        # write the table -
        pretty_table(state_table; header=state_table_header_row)
end
 
# ╔═╡ fd339470-ffef-49fa-8636-dce7924e6405
md"""
### Conclusions
 
The conclusion that we have come to is that unfortunately, this process is not profitable. We tried using different sugars, yet none of them produced a fruitful outcome. No matter how many different methods we tried, we still couldn't produce something with a meaningful profit. We calculated PDGN to be a measley $1.66 per hour while it cost us $552.35 per hour just to run a single reactor. So the cost for running 11 reactors would be $6075.84 per hour while we would only make $18.26 from making PDGN. So we figured that we would make the bare minimum of PDGN while at the same time producing another material that would make us more money. We decided on NAD since we produced the most of that from the reaction. Unfortunately, trying to separate it out using only 8% of the outflow stream and getting it to 95% purity ended up with us making .53 grams of it per hour. That is equivalently worth $62.01. So by running the reactors to max efficieny to maximize profits while still producing enough PDGN to make the buyer happy, we'd end up making $700.37 for every $6075.84 we spent. That doesn't include the inital cost of the reactors, separation units, or pumps which would add more the cost and do nothing to increase the profits. 
"""
 
# ╔═╡ 2f2713eb-a958-4d1a-a1cc-2723ea13c38c
md"""
### References
We got our costs and prices from https://www.sigmaaldrich.com/US/en
"""
 
# ╔═╡ cef22b5d-be5d-49f2-987f-77cf1b9379b9
html"""
<style>
 
main {
    max-width: 1200px;
    width: 75%;
    margin: auto;
    font-family: "Roboto, monospace";
}
 
a {
    color: blue;
    text-decoration: none;
}
 
.H1 {
    padding: 0px 30px;
}
</style>"""
 
# ╔═╡ 213d4486-584f-11ec-2373-5d05e90dc5f8
html"""
<script>
 
        // initialize -
        var section = 0;
        var subsection = 0;
        var headers = document.querySelectorAll('h3, h4');
        
        // main loop -
        for (var i=0; i < headers.length; i++) {
            
               var header = headers[i];
            var text = header.innerText;
            var original = header.getAttribute("text-original");
            if (original === null) {
                
                       // Save original header text
                header.setAttribute("text-original", text);
            } else {
                
                       // Replace with original text before adding section number
                text = header.getAttribute("text-original");
            }
        
            var numbering = "";
            switch (header.tagName) {
                case 'H3':
                    section += 1;
                    numbering = section + ".";
                    subsection = 0;
                    break;
                case 'H4':
                    subsection += 1;
                    numbering = section + "." + subsection;
                    break;
            }
 
               // update the header text 
                header.innerText = numbering + " " + text;
        };
</script>
"""
 
# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BSON = "fbb218c0-5317-5bc6-957e-2ee96dd4b1f0"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
GLPK = "60bf3e95-4087-53dc-ae20-288a0d20c6a6"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
 
[compat]
BSON = "~0.3.4"
DataFrames = "~1.3.0"
GLPK = "~0.15.2"
Plots = "~1.25.1"
PlutoUI = "~0.7.21"
PrettyTables = "~1.2.3"
"""
 
# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised
 
julia_version = "1.6.4"
manifest_format = "2.0"
 
[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "abb72771fd8895a7ebd83d5632dc4b989b022b5b"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.2"
 
[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"
 
[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
 
[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
 
[[deps.BSON]]
git-tree-sha1 = "ebcd6e22d69f21249b7b8668351ebf42d6dc87a1"
uuid = "fbb218c0-5317-5bc6-957e-2ee96dd4b1f0"
version = "0.3.4"
 
[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
 
[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "365c0ea9a8d256686e97736d6b7fb0c880261a7a"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.2.1"
 
[[deps.BinaryProvider]]
deps = ["Libdl", "Logging", "SHA"]
git-tree-sha1 = "ecdec412a9abc8db54c0efc5548c64dfce072058"
uuid = "b99e7846-7c00-51b0-8f62-c81ae34c0232"
version = "0.5.10"
 
[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"
 
[[deps.CEnum]]
git-tree-sha1 = "215a9aa4a1f23fbd05b92769fdd62559488d70e9"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.1"
 
[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f2202b55d816427cd385a9a4f3ffb226bee80f99"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+0"
 
[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "4c26b4e9e91ca528ea212927326ece5918a04b47"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.2"
 
[[deps.ChangesOfVariables]]
deps = ["LinearAlgebra", "Test"]
git-tree-sha1 = "9a1d594397670492219635b35a3d830b04730d62"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.1"
 
[[deps.CodecBzip2]]
deps = ["Bzip2_jll", "Libdl", "TranscodingStreams"]
git-tree-sha1 = "2e62a725210ce3c3c2e1a3080190e7ca491f18d7"
uuid = "523fee87-0ab8-5b00-afb7-3ecf72e48cfd"
version = "0.7.2"
 
[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"
 
[[deps.ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "a851fec56cb73cfdf43762999ec72eff5b86882a"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.15.0"
 
[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"
 
[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"
 
[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "dce3e3fea680869eaa0b774b2e8343e9ff442313"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.40.0"
 
[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
 
[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"
 
[[deps.Crayons]]
git-tree-sha1 = "3f71217b538d7aaee0b69ab47d9b7724ca8afa0d"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.0.4"
 
[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"
 
[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "2e993336a3f68216be91eb8ee4625ebbaba19147"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.0"
 
[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"
 
[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"
 
[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
 
[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
 
[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
 
[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"
 
[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
 
[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"
 
[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"
 
[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"
 
[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"
 
[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"
 
[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"
 
[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"
 
[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"
 
[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"
 
[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
 
[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "0c603255764a1fa0b61752d2bec14cfbd18f7fe8"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.5+1"
 
[[deps.GLPK]]
deps = ["BinaryProvider", "CEnum", "GLPK_jll", "Libdl", "MathOptInterface"]
git-tree-sha1 = "ab6d06aa06ce3de20a82de5f7373b40796260f72"
uuid = "60bf3e95-4087-53dc-ae20-288a0d20c6a6"
version = "0.15.2"
 
[[deps.GLPK_jll]]
deps = ["Artifacts", "GMP_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "fe68622f32828aa92275895fdb324a85894a5b1b"
uuid = "e8aa6df9-e6ca-548a-97ff-1f85fc5b8b98"
version = "5.0.1+0"
 
[[deps.GMP_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "781609d7-10c4-51f6-84f2-b8444358ff6d"
 
[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "30f2b340c2fff8410d89bfcdc9c0a6dd661ac5f7"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.62.1"
 
[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fd75fa3a2080109a2c0ec9864a6e14c60cca3866"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.62.0+0"
 
[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"
 
[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", 
...
