### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 373eac8c-d51c-4ac7-8731-57662d7fda11
begin
	using Pkg
	Pkg.activate("..")
	using TwoLayerDirectNumericalShenanigans, GibbsSeaWater
	using CairoMakie
	using PlutoUI
end

# ╔═╡ 14c81daa-42d8-11ee-1e5c-63e4bd60e8c2
md"""
# Iso-tracer initial condition to compare to marginally stable initial condition

We need to compare the simulations where cabbeling occurs to an exact same setup (i.e. same density profile) where no cabbeling will occur.
This can be done by:

- a *linear* equation of state; or
- an **iso-tracer** profile.

This notebook finds a similar dnesity profile to what we use for the current *large perturbation* experiment but with **isohaline** salinity or **isothermal** temperature.

It might still be the case that we use a linear EOS but at this stage I prefer this option.
"""

# ╔═╡ e8695f4e-31d1-4f99-b634-78f79b32f359
md"""
## Stable system

One option is to leave the salinity field as it as and use an isothermal temperature initial condition.
"""

# ╔═╡ b7f8c64a-cb1c-4a13-9598-f38be1f3243d
begin
	architecture = CPU()
	diffusivities = (ν = 1e-6, κ = (S = 1e-7, T = 1e-7))

	## Setup the model
	model = DNS(architecture, DOMAIN_EXTENT, HIGH_RESOLUTION, diffusivities;
	            reference_density = REFERENCE_DENSITY)

	## set initial conditions
	T₀ᵘ = -1.5
	S₀ᵘ = 34.551
	stable = StableUpperLayerInitialConditions(S₀ᵘ, T₀ᵘ)
	initial_conditions = TwoLayerInitialConditions(stable)
	profile_function = HyperbolicTangent(INTERFACE_LOCATION, 3500.0)

	## `GaussianProfile`
	tracer_perturbation = SalinityGaussianProfile(INTERFACE_LOCATION,
											INTERFACE_LOCATION / 1.1, 100.0, 10.0)

	## With `RandomPerturbations`
	z = znodes(model.grid, Center(), Center(), Center())
	depth_idx = findfirst(z .> INTERFACE_LOCATION / 1.1)
	initial_noise = SalinityNoise(z[depth_idx], 0.001)
	dns = TwoLayerDNS(model, profile_function, initial_conditions; tracer_perturbation, initial_noise)
	set_two_layer_initial_conditions!(dns)
end

# ╔═╡ b3662810-178e-42ff-968d-b33500592d14
fig_ics_stable = visualise_initial_conditions(dns, 1, 1)

# ╔═╡ f5c5ff99-6125-45f6-87c6-f44baf785f23
fig_density_stable = visualise_initial_density(dns, 1, 1, 0)

# ╔═╡ bc86ddff-b42b-4fed-82ed-80ed331e3e33
begin
	S₀ˡ, T₀ˡ = 34.7, 0.5
	Δσ_stable = gsw_sigma0(S₀ˡ, T₀ˡ) - gsw_sigma0(S₀ᵘ, T₀ᵘ)
	md"""
	This is the setup for the **stable** configuration with a *large salinity perturbation* (scale the `GaussianProfile` by ``10``).
	This gives a density diffference between the upper layer (ignoring the perturbation) of ``\Delta\sigma_{0} = \sigma_{0}^{\mathrm{lower}} - \sigma_{0}^{\mathrm{upper}} = `` $(Δσ_stable).
	"""
end

# ╔═╡ f468fbe1-5294-48c6-bff4-269c1b703528
md"""
## Isothermal

One option is to leave the salinity field as it as and use an isothermal temperature initial condition.
"""

# ╔═╡ 1a244403-6b71-40b4-aac0-98e478da773e
begin
	## Setup the model
	model_isothermal = DNS(architecture, DOMAIN_EXTENT, HIGH_RESOLUTION, diffusivities; reference_density = REFERENCE_DENSITY)

	## set initial conditions
	isothermal = IsothermalUpperLayerInitialConditions(34.666, 0.5)
	initial_conditions_isothermal = TwoLayerInitialConditions(isothermal)
	tracer_perturbation_iso = SalinityGaussianProfile(INTERFACE_LOCATION,
											INTERFACE_LOCATION / 1.1, 100.0, 10.0)

	dns_isothermal = TwoLayerDNS(model_isothermal, profile_function, initial_conditions_isothermal, tracer_perturbation = tracer_perturbation_iso; initial_noise)

	## With `RandomPerturbations`
	set_two_layer_initial_conditions!(dns_isothermal)
	fig_ic_isothermal = visualise_initial_conditions(dns_isothermal, 1, 1)
end

# ╔═╡ ce6e0248-9e5f-427f-8f3f-ee92061c889e
fig_density_isothermal = visualise_initial_density(dns_isothermal, 1, 1, 0)

# ╔═╡ 1d6e4a58-a006-48fd-926c-71ed3532524a
gsw_sigma0(S₀ˡ, 0.5) - gsw_sigma0(34.666, 0.5)

# ╔═╡ a5dc9f9b-a28a-4cba-971b-2d702c4420f2
md"""
With a value of 34.666 for the salinity in the upper layer we get a fairly similr density profile and still quite realistic values. Refine this and I think this is the way to go but isohaline option is also below.
"""

# ╔═╡ 29de4ead-4db0-4a56-8819-3e374367d3dd
begin
	S_range = range(16.324, 16.3241, step = 0.0000000001)
	Δσ_iso = gsw_sigma0.(S_range, T₀ˡ) - gsw_sigma0.(S_range, T₀ᵘ)
	find_S = findfirst(Δσ_iso .< Δσ_stable)
	md"""
	## Isohaline
	
	To find an isohaline profile that matches this need to find a salinity value to give us a change in density as above.
	I have already done some playing around and I know that ``S \approx 16.324`` gives a good estimate, ``\Delta\sigma_{0}^{\mathrm{isohaline}} = `` $(gsw_sigma0(16.324, T₀ˡ) - gsw_sigma0(16.324, T₀ᵘ)).

	This salinity value is refined below.
	We have ``S = `` $(S_range[find_S-1]) gives ``\Delta\sigma_{0}^{\mathrm{isohaline}}`` = $(Δσ_iso[find_S-1]).

	This could be refined further but is likely enough for now.
	"""
end

# ╔═╡ 65eb5760-855e-483c-95b6-7fb961ca4fcc
begin
	## Setup the model
	model_iso = DNS(architecture, DOMAIN_EXTENT, HIGH_RESOLUTION, diffusivities;
	            reference_density = REFERENCE_DENSITY)

	## set initial conditions
	isohaline = IsohalineUpperLayerInitialConditions(S_range[find_S-1], T₀ᵘ)
	initial_conditions_iso = TwoLayerInitialConditions(isohaline)
	temperature_perturbation = TemperatureGaussianProfile(INTERFACE_LOCATION,
											INTERFACE_LOCATION / 1.1, 250.0, 1000.0)

	dns_iso = TwoLayerDNS(model_iso, profile_function, initial_conditions_iso, tracer_perturbation = temperature_perturbation)

	## With `RandomPerturbations`
	set_two_layer_initial_conditions!(dns_iso)
	visualise_initial_conditions(dns_iso, 1, 1)
end

# ╔═╡ a94ffb8f-8ab3-48b3-83d6-72f9e3431e6a
fig_density_isohaline = visualise_initial_density(dns_iso, 1, 1, 0)

# ╔═╡ 2a0f4b9e-4df2-4484-8071-e5da25a5f4a3
md"""
# Comparing the two density plots
## Stable
"""

# ╔═╡ 9a54507f-6d07-4d2c-b95b-0c829dcc4b6e
fig_density_stable

# ╔═╡ 45d314d1-8e88-46bb-b81f-d6f5f695b381
md"""
## Isothermal
"""

# ╔═╡ 29b064ae-ee35-48a8-9b8a-d9438b1ce1f3
fig_density_isothermal

# ╔═╡ 55af9d6a-3896-4e88-9676-9df2265b2f84
md"""
## Isohaline
"""

# ╔═╡ 5189781e-31d2-4d21-8cce-cbea9a725d33
fig_density_isohaline

# ╔═╡ e508a976-6173-46ec-b634-581cd7840507
md"""
Up to the small random noise that is seeded (and even this could made the same) at the `GaussianProfile` mean.
"""

# ╔═╡ d65133b1-a33d-4aea-b137-be7e971f362f
TableOfContents(title = "Isohaline density profile")

# ╔═╡ Cell order:
# ╟─14c81daa-42d8-11ee-1e5c-63e4bd60e8c2
# ╟─373eac8c-d51c-4ac7-8731-57662d7fda11
# ╟─e8695f4e-31d1-4f99-b634-78f79b32f359
# ╟─b7f8c64a-cb1c-4a13-9598-f38be1f3243d
# ╟─b3662810-178e-42ff-968d-b33500592d14
# ╟─f5c5ff99-6125-45f6-87c6-f44baf785f23
# ╟─bc86ddff-b42b-4fed-82ed-80ed331e3e33
# ╟─f468fbe1-5294-48c6-bff4-269c1b703528
# ╟─1a244403-6b71-40b4-aac0-98e478da773e
# ╟─ce6e0248-9e5f-427f-8f3f-ee92061c889e
# ╟─1d6e4a58-a006-48fd-926c-71ed3532524a
# ╟─a5dc9f9b-a28a-4cba-971b-2d702c4420f2
# ╟─29de4ead-4db0-4a56-8819-3e374367d3dd
# ╟─65eb5760-855e-483c-95b6-7fb961ca4fcc
# ╟─a94ffb8f-8ab3-48b3-83d6-72f9e3431e6a
# ╟─2a0f4b9e-4df2-4484-8071-e5da25a5f4a3
# ╟─9a54507f-6d07-4d2c-b95b-0c829dcc4b6e
# ╟─45d314d1-8e88-46bb-b81f-d6f5f695b381
# ╟─29b064ae-ee35-48a8-9b8a-d9438b1ce1f3
# ╟─55af9d6a-3896-4e88-9676-9df2265b2f84
# ╟─5189781e-31d2-4d21-8cce-cbea9a725d33
# ╟─e508a976-6173-46ec-b634-581cd7840507
# ╟─d65133b1-a33d-4aea-b137-be7e971f362f
