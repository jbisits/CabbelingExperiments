### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 373eac8c-d51c-4ac7-8731-57662d7fda11
begin
	using Pkg
	Pkg.activate("..")
	using TwoLayerDirectNumericalShenanigans
	using CairoMakie
	using PlutoUI
end

# ╔═╡ 14c81daa-42d8-11ee-1e5c-63e4bd60e8c2
md"""
# Isohaline initial condition to compare to marginally stable initial condition

We need to compare the simulations where cabbeling occurs to an exact same setup (i.e. same density profile) where no cabbeling will occur.
This can be done by:

- a *linear* equation of state; or
- an **isohaline** profile.

This notebook finds a similar dnesity profile to what we use for the current *large perturbation* experiment but with **isohaline** salinity.

It might still be the case that we use a linear EOS but at this stage I prefer this option.
"""

# ╔═╡ d28668b4-ba87-43cb-9e59-9700e558e98f
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
	salinity_perturbation = GaussianProfile(INTERFACE_LOCATION,
											INTERFACE_LOCATION / 1.1, 100.0, 10.0)

	## With `RandomPerturbations`
	z = znodes(model.grid, Center(), Center(), Center())
	depth_idx = findfirst(z .> INTERFACE_LOCATION / 1.1)
	salinity_noise = RandomPerturbations(z[depth_idx], 0.001)
	set_two_layer_initial_conditions!(model, initial_conditions, profile_function,
	                                  salinity_perturbation, salinity_noise)
end

# ╔═╡ b3662810-178e-42ff-968d-b33500592d14
visualise_initial_conditions(model, 1, 1)

# ╔═╡ f5c5ff99-6125-45f6-87c6-f44baf785f23
visualise_initial_density(model, 1, 1, 0)

# ╔═╡ bc86ddff-b42b-4fed-82ed-80ed331e3e33
begin
	Δσ_stable = gsw_sigma0(S₀ˡ, T₀ˡ) - gsw_sigma0(S₀ᵘ, T₀ᵘ)
	md"""
	This is the setup for the **stable** configuration with a *large salinity perturbation* (scale the `GaussianProfile` by ``10``).
	This gives a density diffference between the upper layer (ignoring the perturbation) of ``\Delta\sigma_{0} = \sigma_{0}^{\mathrm{lower}} - \sigma_{0}^{\mathrm{upper}} = `` $(Δσ_stable).

	To find an isohaline profile that matches this need to find a salinity value to give us a change in density as above.
	I have already done some playing around and I know that ``S \approx 16.324`` gives a good estimate, ``\Delta\sigma_{0}^{\mathrm{isohaline}} = `` $(gsw_sigma0(16.324, T₀ˡ) - gsw_sigma0(16.324, T₀ᵘ)).

	This salinity value is refined below.
	"""
end

# ╔═╡ 29de4ead-4db0-4a56-8819-3e374367d3dd
begin
	S_range = range(16.324, 16.3241, step = 0.0000000001)
	Δσ_iso = gsw_sigma0.(S_range, T₀ˡ) - gsw_sigma0.(S_range, T₀ᵘ)
	find_S = findfirst(Δσ_iso .< Δσ_stable)
	md"""
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

	## With `RandomPerturbations`
	set_two_layer_initial_conditions!(model_iso, initial_conditions_iso,
									profile_function, salinity_perturbation, salinity_noise)
	visualise_initial_conditions(model_iso, 1, 1)
end

# ╔═╡ a94ffb8f-8ab3-48b3-83d6-72f9e3431e6a
visualise_initial_density(model_iso, 1, 1, 0)

# ╔═╡ 2a0f4b9e-4df2-4484-8071-e5da25a5f4a3
md"""
# Comparing the two density plots
"""

# ╔═╡ 9a54507f-6d07-4d2c-b95b-0c829dcc4b6e
visualise_initial_density(model, 1, 1, 0)

# ╔═╡ 5189781e-31d2-4d21-8cce-cbea9a725d33
visualise_initial_density(model_iso, 1, 1, 0)

# ╔═╡ e508a976-6173-46ec-b634-581cd7840507
md"""
Up to the small random noise that is seeded (and even this could made the same) at the `GaussianProfile` mean.
"""

# ╔═╡ d65133b1-a33d-4aea-b137-be7e971f362f
TableOfContents(title = "Isohaline density profile")

# ╔═╡ Cell order:
# ╟─14c81daa-42d8-11ee-1e5c-63e4bd60e8c2
# ╟─373eac8c-d51c-4ac7-8731-57662d7fda11
# ╟─d28668b4-ba87-43cb-9e59-9700e558e98f
# ╟─b3662810-178e-42ff-968d-b33500592d14
# ╟─f5c5ff99-6125-45f6-87c6-f44baf785f23
# ╟─bc86ddff-b42b-4fed-82ed-80ed331e3e33
# ╟─29de4ead-4db0-4a56-8819-3e374367d3dd
# ╟─65eb5760-855e-483c-95b6-7fb961ca4fcc
# ╟─a94ffb8f-8ab3-48b3-83d6-72f9e3431e6a
# ╟─2a0f4b9e-4df2-4484-8071-e5da25a5f4a3
# ╟─9a54507f-6d07-4d2c-b95b-0c829dcc4b6e
# ╟─5189781e-31d2-4d21-8cce-cbea9a725d33
# ╟─e508a976-6173-46ec-b634-581cd7840507
# ╟─d65133b1-a33d-4aea-b137-be7e971f362f
