### A Pluto.jl notebook ###
# v0.19.41

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

The main two experiments I will compare are the isothermal experiment and the cabbeling experiment (this is all I have done for DNS).
"""

# ╔═╡ edd033e9-55ca-4b68-bf03-330baaa35e67
md"""
## 1D model - isothermal

This isothermal run of the 1D model has the same setup as the full DNS experiment.
"""

# ╔═╡ d2a1ec90-16ef-47b2-9adb-8829681ad5d9
begin
	iso_data = "../1DModel/OneDModelOutput_isothermal.jld2"
	S = FieldTimeSeries(iso_data, "S")
	T = FieldTimeSeries(iso_data, "T")
	σ₀ = FieldTimeSeries(iso_data, "σ")
	t = S.times
	Δt = diff(t)
	z = znodes(S)
	Δz = S.grid.Δzᵃᵃᶜ
	nothing
end

# ╔═╡ a873d6b3-e6c6-4ed1-ab11-ffaeccb2c4d5
begin
	∫Szdz = sum(interior(S, 1, 1, :, :) .* z * Δz, dims = 1)
	dₜ∫Szdz = vec(diff(∫Szdz, dims = 2))
	fig, ax = lines(t[2:end], dₜ∫Szdz)
	ax.title = "Potential (salt) energy (dₜ∫Szdz)"
	ax.xlabel = "time (s)"
	fig
end

# ╔═╡ 0efe7042-82d9-44aa-be55-eb3634425bd5
begin
	S_sorted = similar(S.data[1, 1, :, :])
	S_content = similar(S.data[1, 1, :, :])
	S_flux = similar(S.data[1, 1, :, 2:end])
	for i ∈ eachindex(t)
		S_sorted[:, i] = sort(interior(S, 1, 1, :, i))	
		S_content[:, i] = cumsum(S_sorted[:, i] * Δz)
	end
	for i ∈ 1:length(t)-1
		S_flux[:, i] = diff(S_content[:, i:i+1], dims = 2) / Δt[i]
	end
end

# ╔═╡ 4a41cc07-e6d2-4aed-b47e-f72e4e03df71
begin
	∫S✶zdz = sum(S_sorted .* z * Δz, dims = 1)
	dₜ∫S✶zdz = vec(diff(∫S✶zdz, dims = 2))
	fig1, ax1 = lines(t[2:end], dₜ∫S✶zdz)
	ax1.title = "Background potential (salt) energy (dₜ∫Szdz)"
	ax1.xlabel = "time (s)"
	fig1
end

# ╔═╡ 944526a6-5154-44ac-af06-c7dba206007c
begin
	fig2, ax2 = lines(t[2:end], dₜ∫Szdz, label = "PE")
	lines!(ax2, t[2:end], dₜ∫S✶zdz, label = "BPE")
	lines!(ax2, t[2:end], dₜ∫Szdz .- dₜ∫S✶zdz, label = "APE")
	ax2.title = "Potential energies"
	axislegend(ax2, position = :rb)
	fig2
end

# ╔═╡ df8c87bc-ba1d-402d-8d40-089ba0f7447f
let
	fig, ax, hm = heatmap(Matrix(S_flux[699:701, :])')
	Colorbar(fig[1, 2], hm)
	fig
end

# ╔═╡ 42f95fcf-be3c-4c8d-8226-5dd33bf73b7a
begin
	S_diff = S_flux[2:end, :] ./ diff(S_sorted[:, 2:end], dims = 1)
	replace!(S_diff, Inf => 0)
	replace!(S_diff, -Inf => 0)
end

# ╔═╡ c5b279df-1031-433f-acf4-9c32c07fcfcb
let
	# fig, ax, hm = heatmap(Matrix(S_diff[695:715, :])')
	# Colorbar(fig[1, 2], hm)
	# fig
	lines(log10.(vec(S_diff[699, :])))
end

# ╔═╡ Cell order:
# ╟─a290de4c-f791-11ee-0e04-11e4e34b7428
# ╟─ac639feb-9ee4-43f4-acd9-466fe3478d40
# ╟─edd033e9-55ca-4b68-bf03-330baaa35e67
# ╟─d2a1ec90-16ef-47b2-9adb-8829681ad5d9
# ╟─a873d6b3-e6c6-4ed1-ab11-ffaeccb2c4d5
# ╟─0efe7042-82d9-44aa-be55-eb3634425bd5
# ╟─4a41cc07-e6d2-4aed-b47e-f72e4e03df71
# ╟─944526a6-5154-44ac-af06-c7dba206007c
# ╠═df8c87bc-ba1d-402d-8d40-089ba0f7447f
# ╠═42f95fcf-be3c-4c8d-8226-5dd33bf73b7a
# ╠═c5b279df-1031-433f-acf4-9c32c07fcfcb
