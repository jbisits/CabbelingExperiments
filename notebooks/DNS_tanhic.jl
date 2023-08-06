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

# ╔═╡ a12a28a5-de8d-4cfb-9ad4-f524cdc002ec
begin
	using Pkg
	Pkg.activate("..")
	using DirectNumericalCabbelingShenanigans
	using DirectNumericalCabbelingShenanigans.TwoLayerDNS
	using DirectNumericalCabbelingShenanigans.OutputUtilities
	using PlutoUI
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
	set_two_layer_initial_conditions!(model, initial_conditions, INTERFACE_LOCATION, :tanh, interface_width; salinity_perturbation = true, salinity_perturbation_width = 100)
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

# ╔═╡ f27bf1d8-407d-44a1-8834-15c868d08b60
visualise_initial_density(model, 0)

# ╔═╡ 2a7a4a77-f597-4f87-82ca-ef3650ca4ffc
md"""
# Look at what random noise looks like
"""

# ╔═╡ 45390e18-6a04-46ff-a2fd-20f2d1176177
add_velocity_random_noise!(model, 1e-20, -0.375)

# ╔═╡ 81f1cbb8-f365-4e57-acab-258214dc9213
model.velocities.u.data[:, :, 4000]

# ╔═╡ 4d54c569-c3d2-410a-b67e-d56ca66b5c1d
findall(znodes(model.grid, Face()) .== -0.375)

# ╔═╡ 512ecbf5-0c3d-4385-9d07-ca732709245c
model.velocities.w.data[:, :, 9]

# ╔═╡ 10c2f596-2285-41f9-b5d4-4d3dc9d78df6
md"""
# Saved data size

Trying to figure out a good save schedule for the data.
Really depends how long the simulation runs for...
"""

# ╔═╡ 6afda154-5286-47c3-860a-ed920bd4abbc
@bind saved_timesteps PlutoUI.Slider(1:2000)

# ╔═╡ 1533101a-47b0-473c-8504-3357404dcca8
begin
	each_save = 20 * 20 * 4000 / 1e9
	saved_bytes = round((each_save) * saved_timesteps; digits = 2)
	text = "Each save step for a given tracer and velocity field is $each_save GB. \nSaving every second for $(round(saved_timesteps / 60; digits = 2)) minutes is $(saved_bytes)GB per tracer and velocity field. \nSo I think that would be around $(saved_bytes)GB * 5 = $(round(saved_bytes * 5; digits = 2))GB."
	Base.Text(text)
end

# ╔═╡ Cell order:
# ╟─a12a28a5-de8d-4cfb-9ad4-f524cdc002ec
# ╟─764fee0a-31bb-11ee-1ec7-8b7bc121a8ed
# ╠═691a96e7-7818-4cd4-8caf-1f70b24e915c
# ╠═86cdf29f-4ca0-4ac2-a673-34d4b8447eec
# ╠═67ac5731-a756-4e5e-9233-eaa9dc36055f
# ╠═f27bf1d8-407d-44a1-8834-15c868d08b60
# ╟─2a7a4a77-f597-4f87-82ca-ef3650ca4ffc
# ╠═45390e18-6a04-46ff-a2fd-20f2d1176177
# ╠═81f1cbb8-f365-4e57-acab-258214dc9213
# ╠═4d54c569-c3d2-410a-b67e-d56ca66b5c1d
# ╠═512ecbf5-0c3d-4385-9d07-ca732709245c
# ╟─10c2f596-2285-41f9-b5d4-4d3dc9d78df6
# ╟─6afda154-5286-47c3-860a-ed920bd4abbc
# ╟─1533101a-47b0-473c-8504-3357404dcca8
