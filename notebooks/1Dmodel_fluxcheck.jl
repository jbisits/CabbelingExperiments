### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ a290de4c-f791-11ee-0e04-11e4e34b7428
begin
	using Pkg
	Pkg.activate("..")
	using JLD2, Oceananigans, CairoMakie
end

# ╔═╡ ac639feb-9ee4-43f4-acd9-466fe3478d40
md"""
# Flux and background state check with 1D model

While computing the background state and the buoyancy flux for the model there were some discrepancies.
In this notebook I run equivalent 1D model and analysis to try and find out what might be going on with some of these calculations.
"""

# ╔═╡ Cell order:
# ╟─a290de4c-f791-11ee-0e04-11e4e34b7428
# ╟─ac639feb-9ee4-43f4-acd9-466fe3478d40
