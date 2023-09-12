### A Pluto.jl notebook ###
# v0.19.27

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

# ╔═╡ 79ad9611-d4b3-4a42-aae7-fd9409f1846c
begin
	using Pkg
	Pkg.activate("..")
	using TwoLayerDirectNumericalShenanigans, CairoMakie, PlutoUI, GibbsSeaWater
end

# ╔═╡ 8f1a175c-503e-11ee-0938-adbbfd497971
md"""
# Initial density profile
For the cabbeling experiments the ideal initial conditions would be one that is outside of the wedge, stable to cabbeling, and one that is *in the wedge, unstable to cabbeling*.
Then just add random noise to kick off the mixing and see what happens.

Initially when I tried to set this as soon as the initial temperature and salinity were such that the profile was unstable to cabbeling the density profile had a bulge.
What we sould like is to be unstable to cabbeling but have a *stable* density profile.

This notebook is exploring how this could be done.
"""

# ╔═╡ 01771446-d79c-4245-8088-52dc5614137e
begin
	architecture = CPU()
	diffusivities = (ν = 1e-6, κ = (S = 1e-7, T = 1e-7))
	
	## Setup the model
	model = DNS(architecture, DOMAIN_EXTENT, HIGH_RESOLUTION, diffusivities;
	            reference_density = REFERENCE_DENSITY)
	model_mp = DNS(architecture, DOMAIN_EXTENT, HIGH_RESOLUTION, diffusivities;
	            reference_density = REFERENCE_DENSITY)
end

# ╔═╡ a1101bd0-55e9-48d3-bcfd-736e1d961bfb
md"""
## Setting salinity and temperature with `HyperbolicTangent`

To determine the stable and unstable to cabbeling initial conditions I have been using my results from project one.
Choosing the stable vs cabbeling initial the instability in the density profile can be seen.

Though this can be worked around by being **less far into the wedge**.
The cabbeling initial condition was halfway between the linearised density and the isopycnal.
If we make this more stable but still in the wedge we can achieve a stable density profile, set by ``S`` and ``T``.
"""

# ╔═╡ 642dbd18-d0aa-46de-9c30-9a3cbb994560
@bind ulic Select(["stable", "cabbeling", "cabbeling2"])

# ╔═╡ dbf28579-232e-440c-815d-6f11607c334c
begin
	T₀ᵘ = -1.5
	S₀ᵘ = (stable = 34.551, cabbeling = 34.568, cabbeling2 = 34.559)
	upper = if isequal(ulic, "stable")
				StableUpperLayerInitialConditions(S₀ᵘ.stable, T₀ᵘ)
			elseif isequal(ulic, "cabbeling")
				CabbelingUpperLayerInitialConditions(S₀ᵘ.cabbeling, T₀ᵘ)
			elseif isequal(ulic, "cabbeling2")
				CabbelingUpperLayerInitialConditions(S₀ᵘ.cabbeling2, T₀ᵘ)
			end
	initial_conditions = TwoLayerInitialConditions(upper)
	profile_function = HyperbolicTangent(INTERFACE_LOCATION, 3500.0)
	depth = find_depth(model, INTERFACE_LOCATION)
	initial_noise = SalinityNoise(depth, 0.001)
	dns = TwoLayerDNS(model, profile_function, initial_conditions; initial_noise)
	set_two_layer_initial_conditions!(dns)
	visualise_initial_density(dns, 1, 1, 0)
end

# ╔═╡ 796798a1-4700-4d1c-9afd-6053a77c7172
let
	S_star, Θ_star = 34.7, 0.5
	S_s, Θ_s = S₀ᵘ.stable, T₀ᵘ
	S_c = isequal(ulic, "cabbeling") ? S₀ᵘ.cabbeling : S₀ᵘ.cabbeling2
	Θ_c = T₀ᵘ
	N = 2000
	S_range, Θ_range = range(34.52, 34.72, length = N), range(-2, 1, length = N)
	S_grid, Θ_grid = ones(N) .* S_range', ones(N)' .* Θ_range
	ρ = gsw_rho.(S_grid, Θ_grid, 0)
	ρ_star = gsw_rho(S_star, Θ_star, 0)
	α_star = gsw_alpha(S_star, Θ_star, 0)
	β_star = gsw_beta(S_star, Θ_star, 0)
	ρ_s = gsw_rho(S_s, Θ_s, 0)
	find_Θ = findfirst(Θ_range .> -1.5)
	find_S = findfirst(ρ[find_Θ, :] .> ρ_star)
	S_iso, Θ_iso = S_range[find_S], Θ_range[find_Θ]
	gsw_rho(S_iso, Θ_iso, 0)

	αₗ, βₗ = gsw_alpha(S_star, Θ_star, 0), gsw_beta(S_star, Θ_star, 0)
	m = βₗ / αₗ
	Θ_linear = @. Θ_star + m * (S_range - S_star)
	fig = Figure(size = (500, 500), fontsize = 22)
	ax = Axis(fig[1, 1];
			  title = "Water masses in salinity-temperature space",
			  xlabel = "Absolute salinity (gkg⁻¹)",
			  ylabel = "Conservative temperature (°C)",
			  limits = (extrema(S_range), extrema(Θ_range)))
	contour!(ax, S_range, Θ_range, ρ'; levels = [ρ_star], color = :black, linewidth = 0.4, labelsize = 18)
	lines!(ax, S_range, Θ_linear, color = :blue, linestyle = :dash)
	scatter!(ax, [S_star], [Θ_star]; color = :red, label = "Lower layer")
	scatter!(ax, [S_s], [Θ_s]; color = :blue, label = "Stable")
	scatter!(ax, [S_c], [Θ_c]; color = :steelblue, label = "Unstable to cabbeling")
	lines!(ax, [S_c, S_star], [Θ_c, Θ_star]; color = :purple, label = "Mixing", linestyle = :dot)
	axislegend(ax, position = :rb)
	fig
end

# ╔═╡ b9910098-4ec6-4e2d-91a4-c8538a91a05c
md"""
In ``S-T`` space this clearly gives us an initial condition that is unstable to cabbeling (the `cabbeling2` initial condition that is) but it is not that unstable to cabbeling. This is a good first test case though and given there is still a fair bit of time on NCI for this quarter two good experiments would be these with only salinity noise to kick off mixing and see what happens.

### Initial S-T profiles that are set
From the figure below I can see that the key thing for setting these initial conditions is **where the S and T nodes are**. If there was a way that the spacing between the salinity and temperature nodes then it could be that the nodes could always be put where they were stable rather than the line joining and having some in the unstable part of ``S-T`` space.
"""

# ╔═╡ 3d801dd0-c1d7-4c6f-9cd1-61cbbbcfdeeb
let
	S_star, Θ_star = 34.7, 0.5
	S_s, Θ_s = S₀ᵘ.stable, T₀ᵘ
	S_c = isequal(ulic, "cabbeling") ? S₀ᵘ.cabbeling : S₀ᵘ.cabbeling2
	S_mid = 0.5*(S_c + S_star)
	Θ_mid = 0.5*(Θ_star + Θ_s)
	T_profile = interior(dns.model.tracers.T, 1, 1, :)
	S_profile = interior(dns.model.tracers.S, 1, 1, :)
	Θ_c = T₀ᵘ
	N = 2000
	S_range, Θ_range = range(34.52, 34.72, length = N), range(-2, 1, length = N)
	S_grid, Θ_grid = ones(N) .* S_range', ones(N)' .* Θ_range
	ρ = gsw_rho.(S_grid, Θ_grid, 0)
	ρ_star = gsw_rho(S_star, Θ_star, 0)
	α_star = gsw_alpha(S_star, Θ_star, 0)
	β_star = gsw_beta(S_star, Θ_star, 0)
	ρ_s = gsw_rho(S_s, Θ_s, 0)
	find_Θ = findfirst(Θ_range .> -1.5)
	find_S = findfirst(ρ[find_Θ, :] .> ρ_star)
	S_iso, Θ_iso = S_range[find_S], Θ_range[find_Θ]
	gsw_rho(S_iso, Θ_iso, 0)

	αₗ, βₗ = gsw_alpha(S_star, Θ_star, 0), gsw_beta(S_star, Θ_star, 0)
	m = βₗ / αₗ
	Θ_linear = @. Θ_star + m * (S_range - S_star)
	fig = Figure(size = (500, 500), fontsize = 22)
	ax = Axis(fig[1, 1];
			  title = "Water masses in salinity-temperature space",
			  xlabel = "Absolute salinity (gkg⁻¹)",
			  ylabel = "Conservative temperature (°C)",
			  limits = (extrema(S_range), extrema(Θ_range)))
	contour!(ax, S_range, Θ_range, ρ'; levels = [ρ_star], color = :black, linewidth = 0.4, labelsize = 18)
	lines!(ax, S_range, Θ_linear, color = :blue, linestyle = :dash)
	scatterlines!(ax, S_profile, T_profile, color = :purple, label = "Initial profile")
	scatter!(ax, [S_star], [Θ_star]; color = :red, label = "Lower layer")
	#scatter!(ax, [S_mid], [Θ_mid])
	axislegend(ax, position = :rb)
	fig
	#save("initialTSprofiles.png", fig)
end

# ╔═╡ 924b4be2-5ed1-4dc1-8f5e-95d5dcbea781
md"""
### Setting a midpoint transition between the upper and lower layers

Can see below that we should solve the problem using a midpoint transition as in figure below.
Clearly this midpoint will only work for *certain* cases but if you can specify where this transition takes pace or find the ``S`` and ``T`` values that give you the *last statically stable transition point* this should work.

As we have a stable profile this should work fine but will have to just set one up and try.
"""

# ╔═╡ d08c96fc-054c-4c0a-836d-da71744f4366
let
	S_star, Θ_star = 34.7, 0.5
	S_s, Θ_s = S₀ᵘ.stable, T₀ᵘ
	S_c = isequal(ulic, "stable") ? S₀ᵘ.stable : isequal(ulic, "cabbeling") ? S₀ᵘ.cabbeling : S₀ᵘ.cabbeling2
	S_mid = 0.5*(S_c + S_star)
	Θ_mid = 0.5*(Θ_star + Θ_s)
	Θ_c = T₀ᵘ
	N = 2000
	S_range, Θ_range = range(34.52, 34.72, length = N), range(-2, 1, length = N)
	S_grid, Θ_grid = ones(N) .* S_range', ones(N)' .* Θ_range
	ρ = gsw_rho.(S_grid, Θ_grid, 0)
	ρ_star = gsw_rho(S_star, Θ_star, 0)
	α_star = gsw_alpha(S_star, Θ_star, 0)
	β_star = gsw_beta(S_star, Θ_star, 0)
	ρ_s = gsw_rho(S_s, Θ_s, 0)
	find_Θ = findfirst(Θ_range .> -1.5)
	find_S = findfirst(ρ[find_Θ, :] .> ρ_star)
	S_iso, Θ_iso = S_range[find_S], Θ_range[find_Θ]
	gsw_rho(S_iso, Θ_iso, 0)

	αₗ, βₗ = gsw_alpha(S_star, Θ_star, 0), gsw_beta(S_star, Θ_star, 0)
	m = βₗ / αₗ
	Θ_linear = @. Θ_star + m * (S_range - S_star)
	fig = Figure(size = (500, 500), fontsize = 22)
	ax = Axis(fig[1, 1];
			  title = "Water masses in salinity-temperature space",
			  xlabel = "Absolute salinity (gkg⁻¹)",
			  ylabel = "Conservative temperature (°C)",
			  limits = (extrema(S_range), extrema(Θ_range)))
	contour!(ax, S_range, Θ_range, ρ'; levels = [ρ_star], color = :black, linewidth = 0.4, labelsize = 18)
	lines!(ax, S_range, Θ_linear, color = :blue, linestyle = :dash)
	scatterlines!(ax, [S_star, S_mid, S_c], [Θ_star, Θ_mid, Θ_c], color = :purple, label = "Initial profile")
	scatter!(ax, [S_star], [Θ_star]; color = :red, label = "Lower layer")
	#scatter!(ax, [S_mid], [Θ_mid])
	axislegend(ax, position = :rb)
	fig
	#save("initialTSprofiles.png", fig)
end

# ╔═╡ a9b3ff4f-c546-494f-9c2d-31664a5aeccc
let
	S_star, Θ_star = 34.7, 0.5
	S_s, Θ_s = S₀ᵘ.stable, T₀ᵘ
	S_c = isequal(ulic, "stable") ? S₀ᵘ.stable : isequal(ulic, "cabbeling") ? S₀ᵘ.cabbeling : S₀ᵘ.cabbeling2
	S_mid = 0.5*(S_c + S_star)
	Θ_mid = 0.5*(Θ_star + Θ_s)
	Δσ = gsw_rho(S_star, Θ_star, 0) - gsw_rho(S_mid, Θ_mid, 0)
end

# ╔═╡ 7d960e45-7257-4283-ac60-89aebaf04976
md"""
## Setting with `MidPoint`

I have written an option to set a `Midpoint` transition at some depth (that can be specified) based on the notes above.

It performs everything as it should so provided it does not blow up when the simulation starts (and this could be refined with more strategically placed nodes) I think it could be there.
"""

# ╔═╡ da31a2a7-f171-417e-a1c4-c1e79b497758
begin
	profile_function_mp = MidPoint(depth)
	dns_mp = TwoLayerDNS(model_mp, profile_function_mp, initial_conditions; initial_noise)
	set_two_layer_initial_conditions!(dns_mp)
	visualise_initial_density(dns_mp, 1, 1, 0)
end

# ╔═╡ 92a878fa-ead4-4cc6-8422-771c2fc8d035
let
	S_star, Θ_star = 34.7, 0.5
	S_s, Θ_s = S₀ᵘ.stable, T₀ᵘ
	S_c = isequal(ulic, "cabbeling") ? S₀ᵘ.cabbeling : S₀ᵘ.cabbeling2
	S_mid = 0.5*(S_c + S_star)
	Θ_mid = 0.5*(Θ_star + Θ_s)
	T_profile = interior(dns_mp.model.tracers.T, 1, 1, :)
	S_profile = interior(dns_mp.model.tracers.S, 1, 1, :)
	Θ_c = T₀ᵘ
	N = 2000
	S_range, Θ_range = range(34.52, 34.72, length = N), range(-2, 1, length = N)
	S_grid, Θ_grid = ones(N) .* S_range', ones(N)' .* Θ_range
	ρ = gsw_rho.(S_grid, Θ_grid, 0)
	ρ_star = gsw_rho(S_star, Θ_star, 0)
	α_star = gsw_alpha(S_star, Θ_star, 0)
	β_star = gsw_beta(S_star, Θ_star, 0)
	ρ_s = gsw_rho(S_s, Θ_s, 0)
	find_Θ = findfirst(Θ_range .> -1.5)
	find_S = findfirst(ρ[find_Θ, :] .> ρ_star)
	S_iso, Θ_iso = S_range[find_S], Θ_range[find_Θ]
	gsw_rho(S_iso, Θ_iso, 0)

	αₗ, βₗ = gsw_alpha(S_star, Θ_star, 0), gsw_beta(S_star, Θ_star, 0)
	m = βₗ / αₗ
	Θ_linear = @. Θ_star + m * (S_range - S_star)
	fig = Figure(size = (500, 500), fontsize = 22)
	ax = Axis(fig[1, 1];
			  title = "Water masses in salinity-temperature space",
			  xlabel = "Absolute salinity (gkg⁻¹)",
			  ylabel = "Conservative temperature (°C)",
			  limits = (extrema(S_range), extrema(Θ_range)))
	contour!(ax, S_range, Θ_range, ρ'; levels = [ρ_star], color = :black, linewidth = 0.4, labelsize = 18)
	lines!(ax, S_range, Θ_linear, color = :blue, linestyle = :dash)
	scatterlines!(ax, S_profile, T_profile, color = :purple, label = "Initial profile")
	scatter!(ax, [S_star], [Θ_star]; color = :red, label = "Lower layer")
	#scatter!(ax, [S_mid], [Θ_mid])
	axislegend(ax, position = :rb)
	fig
	#save("initialTSprofiles.png", fig)
end

# ╔═╡ 3e7fedfd-e0bd-4358-aa46-ce68fd90c015
md"""
## Setting initial density
One option is to set the initial density and temperature profiles then try to reverse engineer the salinity profile to get the stable but unstable to cabbeling profile.

Setting the **density** as the initial condition we can see below that there is no change in the density profile, so it should be possible to reverse engineer a salinity initial condition to match the density profile.
"""

# ╔═╡ 96c44a1f-96d4-4844-a363-f47752b9c014
begin
	density_tanh = HyperbolicTangent(INTERFACE_LOCATION, 3500.0)
	z_model = znodes(dns.model.grid, Center())
	σ₀ˡ = gsw_rho(34.7, 0.5, 0)
	σ₀ᵘ = isequal(ulic, "stable") ? gsw_rho(S₀ᵘ.stable, -1.5, 0) : gsw_rho(S₀ᵘ.cabbeling, -1.5, 0)
	Δσ₀ = σ₀ᵘ - σ₀ˡ
	density_profile = @. σ₀ˡ + 0.5 * Δσ₀ * (1 + tanh(density_tanh.interface_transition_width *
                               (z_model - density_tanh.interface_location)))
	fig, ax = lines(density_profile, z_model)
	ax.xaxisposition = :top
	ax.xlabel = "σ₀ (kgm⁻³)"
	ax.ylabel = "z (m)"
	ax.title = "Initial condition $(ulic)"
	fig
end

# ╔═╡ 9c88315a-11e9-421a-9ec2-f6b1398a5bd3
TableOfContents(title = "Initial density profile")

# ╔═╡ Cell order:
# ╟─8f1a175c-503e-11ee-0938-adbbfd497971
# ╟─79ad9611-d4b3-4a42-aae7-fd9409f1846c
# ╟─01771446-d79c-4245-8088-52dc5614137e
# ╟─a1101bd0-55e9-48d3-bcfd-736e1d961bfb
# ╟─642dbd18-d0aa-46de-9c30-9a3cbb994560
# ╟─dbf28579-232e-440c-815d-6f11607c334c
# ╟─796798a1-4700-4d1c-9afd-6053a77c7172
# ╟─b9910098-4ec6-4e2d-91a4-c8538a91a05c
# ╟─3d801dd0-c1d7-4c6f-9cd1-61cbbbcfdeeb
# ╟─924b4be2-5ed1-4dc1-8f5e-95d5dcbea781
# ╟─d08c96fc-054c-4c0a-836d-da71744f4366
# ╟─a9b3ff4f-c546-494f-9c2d-31664a5aeccc
# ╟─7d960e45-7257-4283-ac60-89aebaf04976
# ╟─da31a2a7-f171-417e-a1c4-c1e79b497758
# ╟─92a878fa-ead4-4cc6-8422-771c2fc8d035
# ╟─3e7fedfd-e0bd-4358-aa46-ce68fd90c015
# ╟─96c44a1f-96d4-4844-a363-f47752b9c014
# ╟─9c88315a-11e9-421a-9ec2-f6b1398a5bd3
