### A Pluto.jl notebook ###
# v0.19.29

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

# ╔═╡ 7c0a9f46-9318-40b7-9765-d5c70b751688
begin
	using Pkg
	Pkg.activate("..")
	using TwoLayerDirectNumericalShenanigans, CairoMakie, PlutoUI, GibbsSeaWater
end

# ╔═╡ 9ac8b892-52a0-11ee-3d54-2771d59c723c
md"""
# `StepChange` profile function

As the `MidPoint` profile function was able to run this notebook is setting up a `StepChange` `profile_function` to check it is all looking good.
Then I want to see if this DNS will run even though there is a discontinuity in the ``S`` and ``T`` profiles that I thought (and might still) cause the simulation to blow up.
"""

# ╔═╡ 613ed422-1bd1-4229-ba32-4ae618e5dd6a
begin
	architecture = CPU()
	diffusivities = (ν = 1e-6, κ = (S = 1e-7, T = 1e-7))
	resolution = (Nx = 10, Ny = 10, Nz = 100)
	## Setup the model
	model = DNS(architecture, DOMAIN_EXTENT, HIGH_RESOLUTION, diffusivities;
	            reference_density = REFERENCE_DENSITY)
end

# ╔═╡ 002c37b8-90b9-4922-b7a6-dac85eb869a6
@bind ulic Select(["stable", "cabbeling"])

# ╔═╡ 6774341e-21f9-40b5-aabb-60593225013d
S₀ᵘ = (stable = 34.551, cabbeling = 34.58)

# ╔═╡ edea3bef-54dc-4dcf-8818-e8691c350888
md"""
## Initial salinity and temperature profiles
"""

# ╔═╡ 55400c9e-63ed-4856-aec6-aa4034514416
begin
	T₀ᵘ = -1.5
	upper = if isequal(ulic, "stable")
				StableUpperLayerInitialConditions(S₀ᵘ.stable, T₀ᵘ)
			elseif isequal(ulic, "cabbeling")
				CabbelingUpperLayerInitialConditions(S₀ᵘ.cabbeling, T₀ᵘ)
			end
	initial_conditions = TwoLayerInitialConditions(upper)
	profile_function = StepChange(INTERFACE_LOCATION)
	depth = find_depth(model, INTERFACE_LOCATION)
	depths = find_depth(model, [INTERFACE_LOCATION + 0.02, INTERFACE_LOCATION - 0.02])
	scales = similar(depths)
	fill!(scales, 2e-4)
	initial_noise = SalinityNoise(depths, scales)
	dns = TwoLayerDNS(model, profile_function, initial_conditions; initial_noise)
	set_two_layer_initial_conditions!(dns)
	visualise_initial_conditions(dns, 1, 1)
end

# ╔═╡ b59482c8-b66d-4182-bb1e-d47829619472
md"""
## Initial density
"""

# ╔═╡ f4ef6c9d-397a-4159-8392-c6db1222d777
visualise_initial_density(dns, 1, 1, 0)

# ╔═╡ 29afdbd7-a662-4640-bcfc-47aef90370b3
md"""
## Initial conditions in ``S-T`` space
"""

# ╔═╡ 1e4eca15-65f3-41f6-870a-5907f09e66f0
let
	S_star, Θ_star = 34.7, 0.5
	S_s, Θ_s = S₀ᵘ.stable, T₀ᵘ
	S_c = S₀ᵘ.cabbeling
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

# ╔═╡ 6bd6898e-f92c-4073-98d1-c14d8d0ca6e0
md"""
## Maximum density

Given that the mixed water *must* fall along the mixing line we should be able to predict the density of the densest water.
I was going to have a look at this theoretically but instead will just use a numerical approach to estimate.
"""

# ╔═╡ 216ddfc4-771a-401d-a805-0395fbe0bc3e
begin
	S_star, Θ_star = 34.7, 0.5
	Θᵘ = -1.5
	slope = (Θᵘ - Θ_star) / (upper.S₀ᵘ - S_star)
	S_mix = range(upper.S₀ᵘ, S_star, step = 0.000001)
	Θ_mix = @. Θᵘ - (slope) * (upper.S₀ᵘ - S_mix)
	#lines(S_mix, Θ_mix)
	ρ_mix = gsw_rho.(S_mix, Θ_mix, 0)
	max_rho, max_rho_idx = findmax(ρ_mix)
	S_max, Θ_max = S_mix[max_rho_idx], Θ_mix[max_rho_idx]
	Δρ_mix = max_rho - gsw_rho(S_star, Θ_star, 0)
	md"""
	Find that the gain in density is $Δρ_mix with the maximum density being $max_rho.
	This is at salinity $(round(S_mix[max_rho_idx], digits = 3))gkg⁻¹ and temperature $(round(Θ_mix[max_rho_idx], digits = 2))°C. So not at midpoint.
	"""
end 

# ╔═╡ 14de00c7-1d40-474f-bbb4-10d3d1baf4ac
let
	S_star, Θ_star = 34.7, 0.5
	S_s, Θ_s = S₀ᵘ.stable, T₀ᵘ
	S_c = S₀ᵘ.cabbeling
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
	ρ_s = gsw_rho(S_s, Θ_s, 0)
	find_Θ = findfirst(Θ_range .> -1.5)
	find_S = findfirst(ρ[find_Θ, :] .> ρ_star)
	S_iso, Θ_iso = S_range[find_S], Θ_range[find_Θ]
	gsw_rho(S_iso, Θ_iso, 0)

	αₗ, βₗ = gsw_alpha(S_star, Θ_star, 0), gsw_beta(S_star, Θ_star, 0)
	m_initial = βₗ / αₗ
	Θ_linear_initial = @. Θ_star + m_initial * (S_range - S_star)
	αₘ, βₘ = gsw_alpha(S_max, Θ_max, 0), gsw_beta(S_max, Θ_max, 0)
	m = βₘ / αₘ
	Θ_linear = @. Θ_max + m * (S_range - S_max)
	fig = Figure(size = (500, 500), fontsize = 22)
	ax = Axis(fig[1, 1];
			  title = "Water masses in salinity-temperature space",
			  xlabel = "Absolute salinity (gkg⁻¹)",
			  ylabel = "Conservative temperature (°C)",
			  limits = (extrema(S_range), extrema(Θ_range)))
	contour!(ax, S_range, Θ_range, ρ'; levels = [ρ_star], color = :red, linewidth = 0.6, labelsize = 18, linestyle = :dot)
	lines!(ax, S_range, Θ_linear_initial, color = :red, linestyle = :dash, label = "Linear density at initial deep water", linewidth = 0.6)
	scatter!(ax, S_profile[1], T_profile[1], color = :red, label = "Deep water")
	scatter!(ax, S_profile[end], T_profile[end]; color = :blue, label = "Shallow water")
	lines!(ax, S_mix, Θ_mix, color = :purple, linestyle = :dot, label = "Mixed water")
	scatter!(ax, S_mix[max_rho_idx], Θ_mix[max_rho_idx], color = :green, label = "Maximum ρ from mixing")
	contour!(ax, S_range, Θ_range, ρ', levels = [max_rho], color = :green, linestyle = :dot, linewidth = 0.6)
	lines!(ax, S_range, Θ_linear, color = :green, linewidth = 0.6, linestyle = :dash, label = "Linear density at new deep water")
	axislegend(ax, position = :lt)
	fig
	#save("initialTSprofiles.png", fig)
end

# ╔═╡ 2dd60351-f416-4beb-97d2-03ca8cb449df
md"""
# Isothermal `StepChange`

As well as running the experiment with the `stable` and `cabbeling` initial conditions above, we will run in the isothermal configuration I used preivously.
This time we will be comparing the `cabbeling` initial condition to an isothermal configuration so the density differnece will be different from last time.

With the `cabbeling` initial condition being more extreme this time the isothermal salinity required to match the initial density difference is *very close* to the salinity in the lower layer.
This is means it is very nearly isothermal and isohaline.
It could mean that the 
"""

# ╔═╡ 153fe311-94f4-4941-94a0-81aecbaec29d
model_isothermal = DNS(architecture, DOMAIN_EXTENT, HIGH_RESOLUTION, diffusivities;
	            reference_density = REFERENCE_DENSITY)

# ╔═╡ 3f0730f5-2f23-4952-8a3a-ed4d1b35ae98
begin
	Δσ_c = gsw_rho(34.7, 0.5, 0) - gsw_rho(34.58, -1.5, 0)
	S_range = range(34.6943, 34.69432, step = 1e-8)
	Δσ_range = gsw_rho(34.7, 0.5, 0) .- gsw_rho.(S_range, 0.5, 0)
	find_s = findfirst(Δσ_range .≤ Δσ_c)
	Δσ_iso = gsw_rho(34.7, 0.5, 0) - gsw_rho(S_range[find_s], 0.5, 0)
	md"""
	If ``S =`` $(S_range[find_s]) then ``\Delta\sigma_{0}^{\mathrm{cab}} ≈ \Delta\sigma_{0}^{\mathrm{isothermal}}`` with $Δσ_c and $Δσ_iso respectively.
	If the accuracy of this needs to be improved can do that.

	## Initial salinity and temperature profiles
	"""
end

# ╔═╡ df8de0eb-57f5-47f4-a330-a2385feb5f46
begin
	upper_isothermal = IsothermalUpperLayerInitialConditions(S_range[find_s], 0.5)
	initial_conditions_isothermal = TwoLayerInitialConditions(upper_isothermal)
	dns_isothermal = TwoLayerDNS(model_isothermal, profile_function, initial_conditions_isothermal; initial_noise)
	set_two_layer_initial_conditions!(dns_isothermal)
	visualise_initial_conditions(dns_isothermal, 1, 1)
end

# ╔═╡ 62209d01-888f-4c5a-adeb-5725691e4993
md"""
## Initial density
"""

# ╔═╡ 66b71108-3210-42bc-ae51-b5abbbb800f4
visualise_initial_density(dns_isothermal, 1, 1, 0)

# ╔═╡ 910bc50d-d8b0-4040-a19c-f0e011480884
md"""
## Initial conditions in ``S-T`` space
"""

# ╔═╡ d958d2ba-e907-4976-9402-64694e331225
let
	S_star, Θ_star = 34.7, 0.5
	S_s, Θ_s = S₀ᵘ.stable, T₀ᵘ
	S_c = S₀ᵘ.cabbeling
	T_profile = interior(dns_isothermal.model.tracers.T, 1, 1, :)
	S_profile = interior(dns_isothermal.model.tracers.S, 1, 1, :)
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

# ╔═╡ 6e41ac9a-84b1-4dfe-bc5c-b3d18b00bac5
TableOfContents()

# ╔═╡ Cell order:
# ╟─9ac8b892-52a0-11ee-3d54-2771d59c723c
# ╟─7c0a9f46-9318-40b7-9765-d5c70b751688
# ╟─613ed422-1bd1-4229-ba32-4ae618e5dd6a
# ╟─002c37b8-90b9-4922-b7a6-dac85eb869a6
# ╟─6774341e-21f9-40b5-aabb-60593225013d
# ╟─edea3bef-54dc-4dcf-8818-e8691c350888
# ╟─55400c9e-63ed-4856-aec6-aa4034514416
# ╟─b59482c8-b66d-4182-bb1e-d47829619472
# ╟─f4ef6c9d-397a-4159-8392-c6db1222d777
# ╟─29afdbd7-a662-4640-bcfc-47aef90370b3
# ╟─1e4eca15-65f3-41f6-870a-5907f09e66f0
# ╟─6bd6898e-f92c-4073-98d1-c14d8d0ca6e0
# ╟─216ddfc4-771a-401d-a805-0395fbe0bc3e
# ╟─14de00c7-1d40-474f-bbb4-10d3d1baf4ac
# ╟─2dd60351-f416-4beb-97d2-03ca8cb449df
# ╟─153fe311-94f4-4941-94a0-81aecbaec29d
# ╟─3f0730f5-2f23-4952-8a3a-ed4d1b35ae98
# ╟─df8de0eb-57f5-47f4-a330-a2385feb5f46
# ╟─62209d01-888f-4c5a-adeb-5725691e4993
# ╟─66b71108-3210-42bc-ae51-b5abbbb800f4
# ╟─910bc50d-d8b0-4040-a19c-f0e011480884
# ╟─d958d2ba-e907-4976-9402-64694e331225
# ╟─6e41ac9a-84b1-4dfe-bc5c-b3d18b00bac5
