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
"""

# ╔═╡ cf0a359e-6b62-4f98-a88c-72bd81e23786
salinity_perturbation = @bind sp Select(["Small", "Medium", "Large"])

# ╔═╡ 1a896dc3-afdf-466f-af1a-865585cbe96d
scale = isequal(sp, "Small") ? 5.0 : isequal(sp, "Medium") ? 8.0 : 10.0

# ╔═╡ f37b8757-2b7c-4539-a30c-10a6013fdab7
begin
	architecture = CPU()
	diffusivities = (ν = 1e-6, κ = (S = 1e-7, T = 1e-7))

	## Setup the model
	model = DNS(architecture, DOMAIN_EXTENT, HIGH_RESOLUTION, diffusivities;
	            reference_density = REFERENCE_DENSITY)
	nothing
end

# ╔═╡ b7f8c64a-cb1c-4a13-9598-f38be1f3243d
begin
	## set initial conditions
	T₀ᵘ = -1.5
	S₀ᵘ = 34.551
	stable = StableUpperLayerInitialConditions(S₀ᵘ, T₀ᵘ)
	initial_conditions = TwoLayerInitialConditions(stable)
	profile_function = HyperbolicTangent(INTERFACE_LOCATION, 3500.0)

	## `GaussianProfile`
	tracer_perturbation = SalinityGaussianProfile(INTERFACE_LOCATION,
											INTERFACE_LOCATION / 1.1, 100.0, scale)

	## With `RandomPerturbations`
	z = znodes(model.grid, Center(), Center(), Center())
	depth_idx = findfirst(z .> INTERFACE_LOCATION / 1.1)
	initial_noise = SalinityNoise(z[depth_idx], 0.001)
	dns = TwoLayerDNS(model, profile_function, initial_conditions; tracer_perturbation, initial_noise)
	set_two_layer_initial_conditions!(dns)
end

# ╔═╡ b3662810-178e-42ff-968d-b33500592d14
fig_ics_stable = visualise_initial_conditions(dns, 1, 1)

# ╔═╡ 0c0b2157-bc39-4ec1-8894-436e1c22c06e
#save("stable_mp_ics.png", fig_ics_stable)

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

# ╔═╡ 852556d9-bba8-4722-a5b3-9c0f27fa688d
begin
	S_isothermalrange = range(34.665, 34.6651, step = 0.0000001)
	Δσ_isothermal = @. gsw_sigma0(S₀ˡ, 0.5) - gsw_sigma0(S_isothermalrange, 0.5)
	find_S_isothermal = findfirst(Δσ_isothermal .< Δσ_stable)
	md"""
	## Isothermal

	One option is to leave the salinity field as it as and use an isothermal temperature initial condition.

	``\Delta\sigma_{0}^{\mathrm{isothermal}} =`` $(Δσ_isothermal[find_S_isothermal]) for salinity value ``S =`` $(S_isothermalrange[find_S_isothermal]).

	With a value of $(S_isothermalrange[find_S_isothermal]) for the salinity in the upper layer we get a fairly similr density profile and still quite realistic values. Refine this and I think this is the way to go but isohaline option is also below.
	"""
end

# ╔═╡ 2d64c18f-f3dc-4859-b9d6-94f90a00ca78
begin
	model_isothermal = DNS(architecture, DOMAIN_EXTENT, HIGH_RESOLUTION, diffusivities; reference_density = REFERENCE_DENSITY)
	nothing
end

# ╔═╡ 1a244403-6b71-40b4-aac0-98e478da773e
begin
	## set initial conditions
	isothermal = IsothermalUpperLayerInitialConditions(S_isothermalrange[find_S_isothermal], 0.5)
	initial_conditions_isothermal = TwoLayerInitialConditions(isothermal)
	tracer_perturbation_iso = SalinityGaussianProfile(INTERFACE_LOCATION,
											INTERFACE_LOCATION / 1.1, 100.0, scale)

	dns_isothermal = TwoLayerDNS(model_isothermal, profile_function, initial_conditions_isothermal, tracer_perturbation = tracer_perturbation_iso; initial_noise)

	## With `RandomPerturbations`
	set_two_layer_initial_conditions!(dns_isothermal)
	fig_ic_isothermal = visualise_initial_conditions(dns_isothermal, 1, 1)
end

# ╔═╡ e6aa5e48-8824-4c48-baf5-9fc66fe044d0
#save("isothermal_mp_ics.png", fig_ic_isothermal)

# ╔═╡ ce6e0248-9e5f-427f-8f3f-ee92061c889e
fig_density_isothermal = visualise_initial_density(dns_isothermal, 1, 1, 0)

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
# ╟─cf0a359e-6b62-4f98-a88c-72bd81e23786
# ╟─1a896dc3-afdf-466f-af1a-865585cbe96d
# ╟─f37b8757-2b7c-4539-a30c-10a6013fdab7
# ╟─b7f8c64a-cb1c-4a13-9598-f38be1f3243d
# ╟─b3662810-178e-42ff-968d-b33500592d14
# ╠═0c0b2157-bc39-4ec1-8894-436e1c22c06e
# ╟─f5c5ff99-6125-45f6-87c6-f44baf785f23
# ╟─bc86ddff-b42b-4fed-82ed-80ed331e3e33
# ╟─852556d9-bba8-4722-a5b3-9c0f27fa688d
# ╟─2d64c18f-f3dc-4859-b9d6-94f90a00ca78
# ╟─1a244403-6b71-40b4-aac0-98e478da773e
# ╠═e6aa5e48-8824-4c48-baf5-9fc66fe044d0
# ╟─ce6e0248-9e5f-427f-8f3f-ee92061c889e
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
