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
    using CairoMakie
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
	profile_function = HyperbolicTangent(INTERFACE_LOCATION, 3500.0)
	set_two_layer_initial_conditions!(model, initial_conditions, profile_function)
	fig = visualise_initial_conditions(model, 1, 1)
end

# ╔═╡ 4829ebc6-ea75-49df-8a95-4ba277f07df1
md"""
## Setting ``\tanh`` steepness

The initial condition for the temperature profiles, similar for salinity, is
```math
\tag{1}
\Theta_{l} + \frac{\Delta \Theta}{2}\left(1 + \tanh s\left(z - i \right)\right)
```
where ``i`` is the `interface_location` and `s` is some scaling to set the steepness of the chnage between the two layers.
What I would like is to be able to **control the number of `znodes` that go over the transition between the two layers.**
If I cannot figure out a function or some other method I can always use this notebook to find the number of nodes that are on the change between the two layers.

The slider `s` below increases ``s`` in ``(1)`` and the plot shows the location of the `znodes`.
In the extreme we see that we get two layers with a step change between them.
In DNS these large jumps (effectively they are discontinuities) will cause the simulation to blow up.
When meeting with Bishakh he said that to ensure the simulation does not blow up require the change to take place over ``10\eta`` where ``\eta`` is the Kolmogorov length scale
```math
\tag{2}
\eta \propto \left(\frac{\nu^{3}}{\epsilon}\right)^{\frac{1}{4}}.
```

I am yet to calculate ``(2)`` ([Oceanostics.jl](https://github.com/tomchor/Oceanostics.jl) can help with ``\epsilon``) but assuming that ``\eta`` is around the size of the vertical resolution in the upper part of the domain want 10`znodes` across the interface.
"""

# ╔═╡ a8cbdf00-42d4-4ab2-8e3d-3d1430cc4887
begin
	s = @bind scale PlutoUI.Slider(1:10000, show_value = true)
	cell_part = @bind cf Select([:Center, :Face])
	zoom = @bind zoomint Select(["Full Profile", "Interface"])
	md"""
	$cell_part
	$zoom
	s = $s
	"""
end

# ╔═╡ 67ac5731-a756-4e5e-9233-eaa9dc36055f
let
	z = cf == :Center ? znodes(model.grid, Center(), Center(), Center()) :
					    znodes(model.grid, Face(), Face(), Face())
	tz = @. 0.5 + (-2 / 2) * (1 + tanh(scale * (z + 0.375)))
	fig, ax = scatterlines(tz, z; label = "z nodes", markersize = 8)
	ax.title = "Location of znodes at $cf for initial condition"
	if isequal(zoomint, "Interface")
		ylims!(ax, -0.425, -0.325)
	end
	axislegend(ax)
	fig
end

# ╔═╡ f27bf1d8-407d-44a1-8834-15c868d08b60
visualise_initial_density(model, 1, 1, 0)

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

# ╔═╡ 01077793-82d9-4222-98a1-590b2b940a96
TableOfContents(title = "Initial conditions")

# ╔═╡ Cell order:
# ╟─a12a28a5-de8d-4cfb-9ad4-f524cdc002ec
# ╟─764fee0a-31bb-11ee-1ec7-8b7bc121a8ed
# ╟─691a96e7-7818-4cd4-8caf-1f70b24e915c
# ╟─86cdf29f-4ca0-4ac2-a673-34d4b8447eec
# ╟─4829ebc6-ea75-49df-8a95-4ba277f07df1
# ╟─a8cbdf00-42d4-4ab2-8e3d-3d1430cc4887
# ╟─67ac5731-a756-4e5e-9233-eaa9dc36055f
# ╟─f27bf1d8-407d-44a1-8834-15c868d08b60
# ╟─10c2f596-2285-41f9-b5d4-4d3dc9d78df6
# ╟─6afda154-5286-47c3-860a-ed920bd4abbc
# ╟─1533101a-47b0-473c-8504-3357404dcca8
# ╟─01077793-82d9-4222-98a1-590b2b940a96
