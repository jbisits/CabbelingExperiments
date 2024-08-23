### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ d7b06e9e-4bf1-11ef-370f-239e7aaa60b0
begin
	using Pkg
	Pkg.activate("..")
	using CairoMakie, JLD2
end

# ╔═╡ fb119988-88f8-46b9-a7d7-e6066f539431
begin
	data = load("../effective_diffusivity.jld2")
end

# ╔═╡ 3ec8d1a5-0134-4ae0-b751-a169b84a802b
let
	fig, ax = lines(data["time/salinity_gradient"][2:end], data["isothermal/salinity_gradient"], label = "isothermal salinity gradient")
	lines!(ax, data["time/salinity_gradient"][2:end], data["cabbeling/salinity_gradient"], label = "cabbeling salinity gradient")
	lines!(ax, data["time/no_gradient"][2:end], data["isothermal/no_gradient"], label = "isothermal no salinity gradient")
	lines!(ax, data["time/no_gradient"][2:end], data["cabbeling/no_gradient"], label = "cabbeling salinity gradient")
	ax.title = "Equal Δρ at model interface"
	ax.xlabel = "time (s)"
	ax.ylabel = "effective diffusivity (m2/s)"
	axislegend(ax, position = :rt)
	fig
end

# ╔═╡ 03f080f4-68de-45bf-8134-84f21a66a68c
let
	fig, ax = lines(data["time/no_gradient"][1:240], data["isothermal/equal_ΔS_cabbeling"])
	ax.xlabel = "time (s)"
	ax.ylabel = "Effective diffusivity (m2/s)"
	ax.title = "Equal ΔS at interface of two layers"
	fig
end

# ╔═╡ b61eb4f9-b498-4643-908f-5d47d5a2b7bf
md"""
Ideas:

- Look at different levels;
- increase salinity gradient;
- save ∂S/∂ẑ
"""

# ╔═╡ 1cae922d-f937-4136-a2a9-4f9fceaa26b3


# ╔═╡ Cell order:
# ╠═d7b06e9e-4bf1-11ef-370f-239e7aaa60b0
# ╠═fb119988-88f8-46b9-a7d7-e6066f539431
# ╟─3ec8d1a5-0134-4ae0-b751-a169b84a802b
# ╟─03f080f4-68de-45bf-8134-84f21a66a68c
# ╠═b61eb4f9-b498-4643-908f-5d47d5a2b7bf
# ╠═1cae922d-f937-4136-a2a9-4f9fceaa26b3
