### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ a12a28a5-de8d-4cfb-9ad4-f524cdc002ec
begin
	using Pkg
	Pkg.activate("..")
	using DirectNumericalCabbelingShenanigans
	using DirectNumericalCabbelingShenanigans.TwoLayerDNS
	using DirectNumericalCabbelingShenanigans.OutputUtilities
end

# ╔═╡ 764fee0a-31bb-11ee-1ec7-8b7bc121a8ed
md"""
# Initial hyperbolic tangent condition for stable DNS profile

This notebook is used to look at what the initially stable profile will look like in the DNS and try to figure out an appropriate thickness for the interface for the experiments (the interface will be the same across all experiments once set).


The model setup is exactly that as in `tanhic_experiments/stable_gpu.jl` which is what will likely be used for the cabbeling experiments.
The main questions that also remain are seeding initial random noise or not and the magnitude of the salinity perturbation in the upper part of the domain.
"""

# ╔═╡ 691a96e7-7818-4cd4-8caf-1f70b24e915c
begin
	architecture = CPU()
	diffusivities = (ν = 1e-4, κ = (S = 1e-5, T = 1e-5))

	## Setup the model
	model = DNS(architecture, DOMAIN_EXTENT, HIGH_RESOLUTION, diffusivities;
	            reference_density = REFERENCE_DENSITY)
end

# ╔═╡ 86cdf29f-4ca0-4ac2-a673-34d4b8447eec
begin
	T₀ᵘ = -1.5
	S₀ᵘ = 34.551
	stable = StableUpperLayerInitialConditions(S₀ᵘ, T₀ᵘ)
	initial_conditions = TwoLayerInitialConditions(stable)
	interface_width = 50
	set_two_layer_initial_conditions!(model, initial_conditions, INTERFACE_LOCATION, :tanh, interface_width; salinity_perturbation = true)
	visualise_initial_conditions(model)
end

# ╔═╡ 67ac5731-a756-4e5e-9233-eaa9dc36055f
let
	T = interior(model.tracers.T, 1, 1, :, 1)
	S = interior(model.tracers.S, 1, 1, :, 1)
	fig, ax = scatterlines(S, T, markersize = 4)
	ax.title = "Initial profile in TS space"
	fig
end

# ╔═╡ 2a7a4a77-f597-4f87-82ca-ef3650ca4ffc
md"""
# Look at what random noise looks like
"""

# ╔═╡ 45390e18-6a04-46ff-a2fd-20f2d1176177
add_veolcity_random_noise!(model, 1e-20, -0.375)

# ╔═╡ 81f1cbb8-f365-4e57-acab-258214dc9213
model.velocities.u.data[:, :, 9]

# ╔═╡ 4d54c569-c3d2-410a-b67e-d56ca66b5c1d
findall(znodes(model.grid, Face()) .== -0.375)

# ╔═╡ 512ecbf5-0c3d-4385-9d07-ca732709245c
model.velocities.w.data[:, :, 9]

# ╔═╡ Cell order:
# ╟─a12a28a5-de8d-4cfb-9ad4-f524cdc002ec
# ╟─764fee0a-31bb-11ee-1ec7-8b7bc121a8ed
# ╠═691a96e7-7818-4cd4-8caf-1f70b24e915c
# ╠═86cdf29f-4ca0-4ac2-a673-34d4b8447eec
# ╠═67ac5731-a756-4e5e-9233-eaa9dc36055f
# ╟─2a7a4a77-f597-4f87-82ca-ef3650ca4ffc
# ╠═45390e18-6a04-46ff-a2fd-20f2d1176177
# ╠═81f1cbb8-f365-4e57-acab-258214dc9213
# ╠═4d54c569-c3d2-410a-b67e-d56ca66b5c1d
# ╠═512ecbf5-0c3d-4385-9d07-ca732709245c
