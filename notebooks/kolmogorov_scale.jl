### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ da38ebd2-601c-4092-bf2f-6ea061e5d0b7
begin
	using Pkg
	Pkg.activate("..")
	using TwoLayerDirectNumericalShenanigans
	using CairoMakie
end

# ╔═╡ fa3ec8c0-471f-11ee-39aa-9b599105dfdd
md"""
# Resolution to satisfy Kolmogorov length scale

From my findings ``\eta \approx 1.6mm``.
The horizontal resolution has been updated to resolve this but need to change grid strecthing and refinement for vertical.

See the plot below but updating `refinement = 1.05` and `stretching = 40` gives a ``\Delta z \approx 7mm`` in the upper 3/4 of the domain.
The domain then increases to a maximum `dz ≈ 1.4mm` which is still less than the Kolmogorov scale and should be ok for out purposes.
"""

# ╔═╡ a306c81f-d6d7-42ff-af76-adda342a683a
begin
	architecture = CPU()
	diffusivities = (ν = 1e-4, κ = (S = 1e-5, T = 1e-5))

	## Setup the model
	model = DNS(architecture, DOMAIN_EXTENT, HIGH_RESOLUTION, diffusivities;
	            reference_density = REFERENCE_DENSITY,
				refinement = 1.05,
				stretching = 40)
end

# ╔═╡ 5847e3cb-9269-43b3-bc9c-1f8773f40459
begin
	zn = znodes(model.grid, Center())
	Δz = zspacings(model.grid, Center(), Center(), Center())
	scatterlines(Δz, zn)
end

# ╔═╡ 0cdf1302-10a5-4873-9aea-ab07bef07d64
Δz[1], Δz[end]

# ╔═╡ 42de73f0-4db0-442d-9247-49812b056d21
findfirst(Δz .== Δz[end])

# ╔═╡ b3c79f92-fa80-430f-ae3b-2eaa366f8708
findfirst(zn .> -0.75)

# ╔═╡ 8256f01f-f9cb-40aa-9e9e-ce4607f00d6b
Δz[300]

# ╔═╡ Cell order:
# ╟─fa3ec8c0-471f-11ee-39aa-9b599105dfdd
# ╟─da38ebd2-601c-4092-bf2f-6ea061e5d0b7
# ╠═a306c81f-d6d7-42ff-af76-adda342a683a
# ╠═5847e3cb-9269-43b3-bc9c-1f8773f40459
# ╠═0cdf1302-10a5-4873-9aea-ab07bef07d64
# ╠═42de73f0-4db0-442d-9247-49812b056d21
# ╠═b3c79f92-fa80-430f-ae3b-2eaa366f8708
# ╠═8256f01f-f9cb-40aa-9e9e-ce4607f00d6b
