### A Pluto.jl notebook ###
# v0.19.42

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
begin
	fig, ax = lines(data["time/salinity_gradient"][2:end], data["isothermal/salinity_gradient"])
	lines!(ax, data["time/salinity_gradient"][2:end], data["cabbeling/salinity_gradient"])
	lines!(ax, data["time/no_gradient"][2:end], data["isothermal/no_gradient"])
	lines!(ax, data["time/no_gradient"][2:end], data["cabbeling/no_gradient"])
	fig
end

# ╔═╡ 03f080f4-68de-45bf-8134-84f21a66a68c
let
	fig, ax = lines(data["time/salinity_gradient"][5:20], data["isothermal/salinity_gradient"][5:20])
	fig
end

# ╔═╡ Cell order:
# ╠═d7b06e9e-4bf1-11ef-370f-239e7aaa60b0
# ╠═fb119988-88f8-46b9-a7d7-e6066f539431
# ╠═3ec8d1a5-0134-4ae0-b751-a169b84a802b
# ╠═03f080f4-68de-45bf-8134-84f21a66a68c
