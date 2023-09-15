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
	depths = find_depth(model, [INTERFACE_LOCATION + 0.01, INTERFACE_LOCATION - 0.01])
	scales = similar(depths)
	fill!(scales, 5e-4)
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

# ╔═╡ 6e41ac9a-84b1-4dfe-bc5c-b3d18b00bac5
TableOfContents()

# ╔═╡ Cell order:
# ╟─9ac8b892-52a0-11ee-3d54-2771d59c723c
# ╟─7c0a9f46-9318-40b7-9765-d5c70b751688
# ╟─613ed422-1bd1-4229-ba32-4ae618e5dd6a
# ╟─002c37b8-90b9-4922-b7a6-dac85eb869a6
# ╠═6774341e-21f9-40b5-aabb-60593225013d
# ╟─edea3bef-54dc-4dcf-8818-e8691c350888
# ╠═55400c9e-63ed-4856-aec6-aa4034514416
# ╟─b59482c8-b66d-4182-bb1e-d47829619472
# ╟─f4ef6c9d-397a-4159-8392-c6db1222d777
# ╟─1e4eca15-65f3-41f6-870a-5907f09e66f0
# ╟─6e41ac9a-84b1-4dfe-bc5c-b3d18b00bac5
